#!/usr/bin/perl -w

# Author: Binary Boy
# This provides my forking system.  Basically, it allows you to queue a
# bunch of tasks.  You set a limit as to the maximum number of tasks to
# run, and it will not allow more than that number to run.

# It runs in the same process and the calling process -- it does not fork
# by itself.  When you create a new instance, it just creates some internal
# data, and nothing more.  When you queue a new command, if this instance
# has not exceeded its limit, then that process is forked and run in the
# background.  Otherwise it is queued (not run).  THE MAIN PROGRAM NEEDS TO
# CALL THE reap() METHOD FROM A REAPER ROUTINE (SIGCHLD handler).  This is
# where the system grabs return values and stuff from the process.  An old
# version of this actually ran more processes from there, but that had
# problems because those processes had SIGCHLD masked (since when the
# fork() occurs, the SIGCHLD is masked since we are still somewhat in the
# signal handler).  In C I can unmask that, but I don't know how to do that
# in Perl, and I'm not sure if it's a good idea anyways.  So, instead,
# several other methods in this module check to see if the tasks are fully
# populated and run more if necessary.  Thus, in order to run fairly
# efficiently, you need to call things from this method occassionally.
# Hopefully that's not too bad -- if you run something, you're probably
# interested in some result of it, so you'll at least want to know when it
# finished, and the method to check that will do that.  I mention which
# ones do that in their descriptions below.  Tasks are run in a FIFO
# manner.  Tasks have their stdin, stdout, and stderr rerouted, so they
# never connect to the terminal in any fashion.  However, you can get the
# output of their stdout and stderr (see the command structure part below).

# It should be possible to have multiple instances of this class, and each
# one will keep track of its processes, making sure not to go over a limit.
# That way, you could have one for controlling downloading from certain
# sites, and another for controlling downloading from other sites, and have
# different limits as to the number of simultaneous downloads from each
# site.

# An interesting race condition that I ran into with the previous version:
# The reap() method, called from a reaper routine, would start new
# processes as space opened up.  You could call it directly to start
# processes, too.  If you did that, if was possible for the new process to
# finish (especially if the given command didn't exist) before the new one
# had recorded the PID, so reap() would get called again by the signal
# handler, and would ignore it since the PID wasn't in its list.  Then, it
# would go back to the original instance of reap(), which would record the
# PID.  That's a PID that would never finish, since it didn't exist
# anymore.  That doesn't happen if you ONLY call reap() from the signal
# handler (since signals are blocked until one is handled), but that's a
# poor solution, since the subsequent processes had the signals blocked,
# too.

# FOR THE REASON DESCRIBED IN THE ABOVE PARAGRAPH, DO NOT RUN METHODS OTHER
# THAN reap() IN YOUR SIGNAL HANDLER.  Other methods may spawn processes,
# and those processes will not be able to recieve SIGCHLD.  If you're OK
# with that, go for it, but you have been warned.


{ package Fork_System;
use strict;
use warnings;
# From the perlipc man page, this is for WNOHANG.
use POSIX ();


# Constructor:
#   new_fork() -- the constructor.  Just creates a new data structure, and
#     returns it.  Both of these should work:
#       Fork_System->new_fork()
#       Fork_System::new_fork().
#     Even this should work:
#       $fork_system->new_fork()
#     All of these should return a new instance (blessed hash ref), of the
#     appropriate class.

# Methods:
#   running() --
#     Returns true if there are processes queued or running, false
#     otherwise.  Attempts to run queued commands if there is room.
#   status() --
#     Returns a two element array, containing the number of queued tasks
#     and the number of running tasks, respectively.  Attempts to run
#     queued commands if there is room.
#   verbosity(verbosity) --
#     Sets and returns the verbosity of this class.  If no argument is
#     given, it sets nothing (just returns).  If an argument is given, the
#     verbosity will be set, and the previous value returned.  verbosity
#     controls how much it dumps out stdout (nothing is ever sent to stderr
#     from here).  The levels are numeric, like so:
#     	0: no output
#     	1: error messages only
#     	2: informational messages
#     	3 and above: various debugging levels
#   max(max) --
#     Sets and returns the maximum number of processes that this instance
#     will run at a time.  With no argument, nothing is changed, and the
#     current value is returned.  With an argument, the value is set, and
#     the previous value returned.  If you increase it, new tasks will not
#     be run until reap() is called or a new task is queued.  To get
#     immediate response, try running reap() with a PID of 0, which will be
#     ignored (or no argument at all).  However, this leads to a race
#     condition (described under the reap() method), so the recommended
#     way of doing this is to call it with: kill("CHLD",$$).  If you
#     decrease the max number, the number of running tasks will NOT be
#     decreased (which would require killing stuff).  Rather, the routines
#     will just let tasks finish, and not start new ones until the number
#     drops below this level.
#   queue(com) --
#     This queues a command.  The com argument is a command structure (see
#     below).  This returns 1 if the process was successfully queued, and
#     undef for some error (not a valid command or something).  To find out
#     if the command is running, check the "pid" element of the command
#     structure (if you defined it).
#   reap(pid,retval) --
#     This needs to be called from the main program in the reaper routine
#     for prompt response.  It is OK to call it with a PID that does not
#     belong to it.  That's quite all right.  This alters the internal
#     data, including storing return values, closing streams, etc. so that
#     the calling process can see what happened just by looking at the
#     command structure.
#   killall(signal) --
#     Sends the specified signal (SIGTERM by default) to all the running
#     children of this fork.  This is probably best used for killing all
#     children when you want to shut down.  Try something like:
#       $fork->max(0);  # Prevent new tasks from starting.
#       my($num_run,$old_num_run);
#       (undef,$old_num_run)=$fork->status();
#       $fork->killall();  # Send kill signal.
#       # Wait a bit -- a simple sleep() will get interrupted as stuff dies.
#       for(;;){
#         sleep(1);
#         (undef,$num_run)=$fork->status();
#         if($num_run>=$old_num_run){last;}
#         $old_num_run=$num_run;
#       }
#       $fork->killall(9);
#     If you need something special done, you had best keep track of PIDs
#     and do the killing yourself.

# Destructor:
#   No destructor.  When this passes out of scope, the children continue to
#   run.  It is probably best to not let it pass out of scope until you are
#   sure all the children are done.  You might encourage that by keeping
#   track of all the processes you started, and then killing them when you
#   want to finish up (or use killall() -- see the comments for that method
#   above).

# All methods set the hash element "error" to an error message in the event
# of an error (but nothing ever clears it).  Feel free to check that in the
# event of an error.

# Usage note:  The calling program should use the methods listed above, and
# not call any other methods.  Also, looking at or changing the internal
# variables void the warranty.  You may look at, but do not change, the
# command structures you pass to the queue routine.

# Command structures:
# These are hash references, which are passed to the queue routine and
# stored internally.  The hash needs to contain a key named "command",
# which is a scalar that is the command to be executed.  command is passed
# to exec directly, and can be a string or an array reference (in which
# case, the array is passed to exec).  See "indirect" below.
# Optional keys are:
#   "stdin", "stdout", "stderr" (which will be filehandles for the
#     associated inputs/outputs -- note that if these are defined, the
#     calling process is responsible for monitoring them, so that the
#     processes don't hang waiting for input, or for processing output.
#     The same goes for closing them when done.  If not given, /dev/null
#     will be used.  If set to undef, they will not be reopened.)
#   "retval" (which will contain the return value of the command, and will
#     be set to -1 until the command is completed)
#   "pid" (which will contain the PID when running, 0 before, and deleted
#     afterwards)
#   "indirect" (this is passed to exec as an indirect object, which can
#     override the command that will actually run, and force it to be not
#     interpreted by the shell).
#   "path" (if this exists, the child will attempt to chdir to "path"
#     before running).
# NOTE: these are stored within the module only as long as the task is
# running.  Once removed, the calling program, if it does not keep the
# reference it passed to queue(), loses all the data in them.  However, I
# remind you again, the parent is not encouraged to alter the data in the
# structure.


# FORKING SYSTEM INTERNAL CODE ########
# These routines are not meant to be called by the main program, and should
# not be exported.


# Given a command structure, this prints out the command part in a
# read-able fashion.  Returns an empty string for success, and an error
# message for failure.
sub print_com($){
	my($comref)=@_;

	if(ref($comref) ne "HASH"){
		return "print_com(): not given a hash reference.";
	}
	if(!defined($comref->{"command"})){
		return "print_com(): no command key.";
	}
	if(defined($comref->{"indirect"})){
		print "{".$comref->{"indirect"}."} "
	}
	if(ref($comref->{"command"}) eq "ARRAY"){
		local $,=" ";print @{$comref->{"command"}};
	}else{print $comref->{"command"};}
	return "";
}


# Given a command structure, this forks a child, sets up handles,
# initializes the data, and then runs the command.  Returns an error string
# (will start with "run_com():") for error, and the PID for success (or, so
# it thinks -- the child may have some issues getting started, and this
# won't notice, but that is not likely).  The optional second argument is a
# reference to the running_task hash.  If it's there, when this forks, it
# will record the command reference under the PID of the task in the hash.
sub run_com($;$){
	my($comref,$runref)=@_;
	my($handle,$stdin,$stdout,$stderr);

	if(ref($comref) ne "HASH"){
		return "run_com(): not given a hash reference.";
	}
	if(!defined($comref->{"command"})){
		return "run_com(): no command key.";
	}

	# Prepare the requested pipes.
	if(defined($comref->{"stdin"})){
		undef($comref->{"stdin"});
		if(!pipe($stdin,$comref->{"stdin"})){
			return "run_com(): stdin pipe: $!";
		}
		binmode($stdin);binmode($comref->{"stdin"});
		$handle=select($stdin);$|=1;select($comref->{"stdin"});$|=1;
		select($handle);
	}
	if(defined($comref->{"stdout"})){
		undef($comref->{"stdout"});
		if(!pipe($comref->{"stdout"},$stdout)){
			if(defined($stdin)){close($stdin);close($comref->{"stdin"});}
			return "run_com(): stdout pipe: $!";
		}
		binmode($stdout);binmode($comref->{"stdout"});
		$handle=select($stdout);$|=1;select($comref->{"stdout"});$|=1;
		select($handle);
	}
	if(defined($comref->{"stderr"})){
		undef($comref->{"stderr"});
		if(!pipe($comref->{"stderr"},$stderr)){
			if(defined($stdin)){close($stdin);close($comref->{"stdin"});}
			if(defined($stdout)){close($stdout);close($comref->{"stdout"});}
			return "run_com(): stderr pipe: $!";
		}
		binmode($stderr);binmode($comref->{"stderr"});
		$handle=select($stderr);$|=1;select($comref->{"stderr"});$|=1;
		select($handle);
	}

	# Fork.
	$handle=fork();if(!defined($handle)){
		return "run_com(): fork() failed: $!";
	}
	if($handle>0){
		# Parent -- record a new running task, set PID if necessary and return
		# success, to the best of our knowledge.
		if(defined($runref) && ref($runref) eq "HASH"){$runref->{$handle}=$comref;}
		if(defined($comref->{"pid"})){$comref->{"pid"}=$handle;}
		return $handle;
	}
	# Child -- attempt to get into a different process group (but stay in
	# same session) so that terminal things like CTRL-C don't make it to
	# here, set handles, and run like mad.
	# An error while reopening STDERR results in a silent death.  Anything
	# else may be safely sent to STDERR.
	setpgrp(0,0);
	if(defined($stderr)){open(STDERR,">&",$stderr) || exit(1);}
	elsif(!exists($comref->{"stderr"})){open(STDERR,">","/dev/null") || exit(1);}
	if(defined($stdin)){open(STDIN,"<&",$stdin) || die "Could not reopen STDIN.\n";}
	elsif(!exists($comref->{"stdin"})){
		open(STDIN,"<","/dev/null") || die "Could not reopen STDIN.\n";
	}
	if(defined($stdout)){open(STDOUT,">&",$stdout) || die "Could not reopen STDOUT.\n";}
	elsif(!exists($comref->{"stdout"})){
		open(STDOUT,">","/dev/null") || die "Could not reopen STDOUT.\n";
	}
	if(defined($comref->{"path"})){
		chdir($comref->{"path"}) || die "Could not chdir to ".$comref->{"path"}.".\n";
	}
	if(defined($comref->{"indirect"})){
		if(ref($comref->{"command"}) eq "ARRAY"){
			( exec {$comref->{"indirect"}} (@{$comref->{"command"}}) ) ||
				die "Could not run exec().\n";
		}else{
			( exec {$comref->{"indirect"}} ($comref->{"command"}) ) ||
				die "Could not run exec().\n";
		}
	}else{
		if(ref($comref->{"command"}) eq "ARRAY"){
			exec(@{$comref->{"command"}}) || die "Could not run exec().\n";
		}else{
			exec($comref->{"command"}) || die "Could not run exec().\n";
		}
	}
}


# This populates the fork system by running queued tasks until the queue is
# empty of the max number of processes are running.
# Needs the hash reference to this instance of the class.
# Returns a true for success, and an error message for failure (will start
# with "populate():").
sub populate($){
	my($this)=@_;

	# This needs to be a valid instance.
	if(ref($this) eq ""){return "populate(): invalid instance.";}

	# Cut down on number of dereferences with these.
	my $running_task=$this->{"running_task"};
	my $task_queue=$this->{"task_queue"};
	my $verbosity=$this->{"verbosity"};
	my $com;

	# Now, run tasks from the queue until it is empty, or until we have as
	# many tasks running as we want.
	# Re-use $pid -- I have no more use for its value.
	while(scalar(keys(%{$running_task}))<$this->{"max"}){
		# Use shift here.  If you push stuff on, and shift if off, you get a
		# FIFO.
		$com=shift(@{$task_queue});if(!defined($com)){
			if($verbosity>=3){print "populate(): no tasks queued -- not running more.\n";}
			last;
		}
		# The hidden assumption here is that I will be able to add the PID to
		# my internal list before I get another SIGCHLD, but I do not know how
		# to guarantee that save that I get to that very quickly.  Also, I hope
		# perl queues signals so that as long as this routine was called within
		# the SIGCHLD handler, then it won't be called again for the next
		# SIGCHLD until this one finishes.  IMPORTANT:  Only call it from the
		# SIGCHLD handler (and signal yourself to force its invocation).
		my $pid=run_com($com,$running_task);if($pid=~/^run_com\(\):/){
			my $msg="populate(): $pid";
			if($this->{"verbosity"}>=1){print $msg."\n";}
			return $msg;
		}
		if($this->{"verbosity"}>=2){
			print "populate(): started $pid.  Command:\n  ";
			print_com($com);print "\n";
		}
	}

	# Done.
	return 1;
}


# FORKING SYSTEM FRONT END CODE #######
# These are the routines that the main program is supposed to use.


# Here is the constructor.
sub new_fork(;$){
	my($class)=@_;

	my %self=();my $this=\%self;

	# Initialize a bunch of default values.
	$self{"max"}=5;$self{"verbosity"}=0;
	$self{"running_task"}={};$self{"task_queue"}=[];

	# Now, return.  If we have a specific class to be a part of (if this was
	# called using the "->" operator), use it ...
	if(defined($class)){
		if(ref($class) ne ""){return bless($this,ref($class));}
		return bless($this,$class);
	}
	# ... otherwise use a default.
	return bless($this);
}


# This removes the given process from the running task list (the assumption
# is that you have a reaper process calling this with the PID of processes
# that have died).  Returns 1 for success, undef for error.
sub reap($;$$){
	my($this,$pid,$retval)=@_;

	# This needs to be a valid instance.
	if(ref($this) eq ""){return undef;}

	if($this->{"verbosity"}>=3){
		print "reap(): called with pid=";
		if(!defined($pid)){print "<none>";}else{print $pid;}
		print " and retval=";
		if(!defined($retval)){print "<none>";}else{print $retval;}
		print ".\n";
	}

	# Cut down on number of dereferences with these.
	my $running_task=$this->{"running_task"};
	my $task_queue=$this->{"task_queue"};
	my $verbosity=$this->{"verbosity"};
	my $com;

	# First, remove from the running task list, if appropriate.
	if(defined($pid) && defined($running_task->{$pid})){
		if($verbosity>=2){print "reap(): removing $pid from running_task.\n";}
		# Go ahead and change the data, because, even though we delete the
		# reference here, there may be a reference to it elsewhere.
		if(!defined($retval)){
			if($verbosity>=3){print "reap(): warning: no retval was passed.\n";}
		}elsif(defined($running_task->{$pid}->{"retval"})){
			$running_task->{$pid}->{"retval"}=$retval;
		}
		delete $running_task->{$pid};
	}

	return 1;
}


# This queues the given task (or runs it, if that's OK).
sub queue($$){
	my($this,$com)=@_;

	# This needs to be a valid instance.
	if(ref($this) eq ""){return undef;}

	if(ref($com) ne "HASH"){
		$this->{"error"}="queue(): Called without a valid command structure.";
		if($this->{"verbosity"}>=1){print $this->{"error"}."\n";}
		return undef;
	}

	if(!defined($com->{"command"})){
		$this->{"error"}="queue(): No command key command structure.";
		if($this->{"verbosity"}>=1){print $this->{"error"}."\n";}
		return undef;
	}

	# Set the stuff.
	if(defined($com->{"pid"})){$com->{"pid"}=0;}
	if(defined($com->{"retval"})){$com->{"retval"}=-1;}
	# Add it to the queue.
	push(@{$this->{"task_queue"}},$com);
	if($this->{"verbosity"}>=2){
		print "Added to command queue:\n  ";
		print_com($com);print "\n";
	}
	if($this->{"verbosity"}>=3){
		my @list=();local $,="";
		if(defined($com->{"indirect"})){print "  Indirect: ",$com->{"indirect"},"\n";}
		if(defined($com->{"path"})){print "  Path: ",$com->{"path"},"\n";}
		if(defined($com->{"stdin"})){push(@list,"stdin");}
		if(defined($com->{"stdout"})){push(@list,"stdout");}
		if(defined($com->{"stderr"})){push(@list,"stderr");}
		if(defined($com->{"retval"})){push(@list,"retval");}
		if(defined($com->{"pid"})){push(@list,"pid");}
		local $,=" ";print "  Other defines:",@list,"\n";
	}
	# Run any commands necessary.
	my $ret=populate($this);if($ret=~/^populate\(\):/){
		$this->{"error"}="queue(): $ret";
		if($this->{"verbosity"}>=1){print $this->{"error"}."\n";}
		return undef;
	}
	return 1;
}


sub running($){
	my($this)=@_;

	# This needs to be a valid instance.
	if(ref($this) eq ""){return undef;}

	# Run any commands necessary.
	my $ret=populate($this);if($ret=~/^populate\(\):/){
		$this->{"error"}="running(): $ret";
		if($this->{"verbosity"}>=1){print $this->{"error"}."\n";}
		return undef;
	}

	if(scalar(@{$this->{"task_queue"}}) ||
	   scalar(keys(%{$this->{"running_task"}}))){return 1;}
	return 0;
}


sub status($){
	my($this)=@_;

	# This needs to be a valid instance.
	if(ref($this) eq ""){return undef;}

	# Run any commands necessary.
	my $ret=populate($this);if($ret=~/^populate\(\):/){
		$this->{"error"}="running(): $ret";
		if($this->{"verbosity"}>=1){print $this->{"error"}."\n";}
		return undef;
	}

	return (scalar(@{$this->{"task_queue"}}),
	        scalar(keys(%{$this->{"running_task"}})));
}


sub verbosity($;$){
	my($this,$v)=@_;

	# This needs to be a valid instance.
	if(ref($this) eq ""){return undef;}

	my $prev=$this->{"verbosity"};
	if(defined($v)){$this->{"verbosity"}=$v;}
	if($this->{"verbosity"}>=2){
		local $,="";
		print "verbosity(): verbosity set to ",$this->{"verbosity"},"\n";
	}
	return $prev;
}


sub max($;$){
	my($this,$m)=@_;

	# This needs to be a valid instance.
	if(ref($this) eq ""){return undef;}

	my $prev=$this->{"max"};
	if(defined($m)){$this->{"max"}=$m;}
	if($this->{"verbosity"}>=2){
		local $,="";
		print "max(): max set to ",$this->{"max"},"\n";
	}
	return $prev;
}


sub killall($;$){
	my($this,$sig)=@_;

	# This needs to be a valid instance.
	if(ref($this) eq ""){return undef;}

	if(!defined($sig)){$sig="TERM";}

	kill($sig,keys(%{$this->{"running_task"}}));
}


# End of Fork_System class.
}

# So that the loading is deemed successful.
return 1;
