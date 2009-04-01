#!/usr/bin/perl -w


# A library of calls for use in our experiments.
# Exports several function calls into the caller namespace.

# NOTE: This library forks(), and sometimes expects to receive SIGCHLD.  It
# traps that signal on initialization.  If you wish to take it over, be
# sure to call the reap() routine within this library with the PID and
# return value of every process that may be related to this library (best
# to err on the side of passing ALL processes) in order for it to function
# properly.  reap() is NOT exported by default.

# NOTE: The above mentioned property has a strange side effect: system()
# returns failure, because it does not catch the child dying (reports error
# "no child processes" or something like that).  To counter this, use
# run_com(), queue_com(), and wait_com() in this module (or just don't care
# about that particular error) (also available with *_silent versions).

# Also, one little thing:  there's a non-exported flag, Experiment::testing.
# Normally, it is set to 0.  If you want to test a programming setup, set
# it to one, and the program_boxes() routine will generate and save
# programs, but will NOT connect to boxes, send them data, UDP ping them,
# or whatever.  Also, the trigger() and wait_for_boxes() routines will just
# return immediately.  Likewise, the cam_takepic() and cam_wait() return
# immediately without doing anything.  You can still make the camera wait
# for an external trigger with the cam_command() routine, so I recommend
# not doing that.  This way, you can talk to the camera even in testing
# mode.

# This library had lot of changes from the previous version -- it allows
# some way to have add-ons, and should now initialize the outputs on analog
# boxes so that the initial conditions don't have to be set at the
# beginning of the programs.

# Also, I fixed a conversion problem.  On the analog boxes, the integer
# code 0xFFFF is NOT 10V, like I thought.  0x10000 is 10V, which is not
# quite achievable, so I modified the conversions to correct for that.


BEGIN {
package Experiment;

use strict;
use warnings;
use Socket;
use Fcntl;
use POSIX;
use Time::HiRes;
use Carp;
use Exporter qw(import);
use Fork_System;


# GLOBALS #############################


# The special testing flag mentioned in the comments in the header:
our $testing=0;


# This library is based on a data type called an Action.  Each type of box
# it controls has its own Action.  Digital boxes have Digital_Action,
# analog boxes have Analog_Action, etc.  I wanted to make it possible to
# add a separate module that added another action, so here's where you add
# it on.  Any add-on module would need to append itself to this hash.  The
# keys are the names of Actions.  You provide code that a user can call to
# produce Actions, which are then stuck together in arrays of actions
# and/or Events.  When the user calls program_boxes(), it references this
# hash for every Action type it finds.  It then does the following (buried
# in here are the various hash elements it needs):
# 1) Creates the data needed for the box from the actions:
#      'prog_data' (required) -- a subroutine reference that creates the
#        data.  See the digital_box_data() or analog_box_data() routines
#        for the arguments, return values, etc.  It also sets a value
#        determining whether the box needs to be UDP pinged.  The one
#        exception is the master box, which is never pinged by
#        program_data() (it may be pinged by trigger(), which runs off a
#        different criteria).
# 2) Initializes the box:
#      'init' (optional) -- a subroutine reference that sets up the box the
#        way it should be before starting the run.  See analog_box_init()
#        for the comments.  If 'init' is defined, it will get called when
#        the time comes to program boxes.  It is the responsibility of the
#        routine to recognize whether or not there's data there, and
#        whether or not it should actually do something.
#      Initializations are performed one box at a time.  This is because
#      the initializations are assumed to be somewhat uncommon (few boxes
#      need them), not too slow, and complicated enough that it would be
#      hard to interrupt and go service something else.
#      Initializations are always performed if the data is there.  If the
#      main programming is sufficiently complicated, you can have it done
#      here as well.
# 3) Programs the box:
#      'port' (optional) -- the TCP port number to which to connect for
#        sending program data.  If this is undefined, or no program data
#        was created by 'prog_data', then no program is sent.
#      'respond' (required if 'port' is defined) -- a routine that sends
#        the program data to the box, responding to its requests.  See
#        respond_to_box() for an example.
#      program_boxes() opens a bunch of sockets at a time, and monitors all
#      of them, calling the appropriate 'respond' every time it receives
#      something from the boxes.  This is because all the programming we do
#      is pretty basic request/response-type stuff that is easily
#      state-machined, and the computer is capable of doing things much
#      faster than the Ethernuts on the other end (as of this writing).
#      Boxes are only programmed if the program data is different than what
#      it was the last time the box was programmed by this module.  If this
#      is the first time programming the box, and this is the first time
#      the module has programmed it (since being loading into memory), then
#      the box will be programmed.
# 4) Sends UDP pings to boxes to prepare them for the trigger.  The master
#    box is not pinged until the trigger() routine is called (outside of
#    program_boxes()).  The 'prog_data' routine should set a 'UDP' element
#    for the box data to some true value if it wants that particular box
#    UDP pinged.
# There is one more thing that is performed when creating Event objects.
# That requires copying data.  Simple copies can be done by the existing
# code, but if you have references to stuff in the Action object, you'll
# get multiple references to the same data.  Then, changing that data will
# change EVERY Event that references the object.  If you don't want that,
# add one more element:
#      'copy' (optional) -- if it exists, then it is a reference to a
#        routine used to copy Actions.  See dig_act_copy() for an example.
# This arrangement makes for pretty messy code within this module, but it
# should be OK if you are just writing stuff outside this module.  It could
# probably have been cleaned up a little by using objects, but that would
# require people adding modules to know a bit more about objects in Perl,
# and it seemed easier to just leave it be.
our %Action=('Digital_Action'=>{'prog_data'=>\&digital_box_data,
                                'respond'  =>\&respond_to_box,
                                'init'     =>undef,
                                'copy'     =>\&dig_act_copy,
                                'port'     =>24,
                               },
             'Analog_Action' =>{'prog_data'=>\&analog_box_data,
                                'respond'  =>\&respond_to_box,
                                'init'     =>\&analog_box_init,
                                'copy'     =>undef,
                                'port'     =>24,
                               },
);


# The list of parameters to record with each run.
our @parameters=();
# Name of parameters file.
our $PARAMETERS='PARAMETERS';


# Each key is the name of a box, and the value, if defined, is the string
# that the box was last programmed with.
our %programs=();
# The most recently programmed boxes:
our @programmed_boxes=();


# When you initialize a run (see init_run()), a directory is made, we chdir
# into it (recording the dir we were in in here), save parameters, etc.
# When we end a run (see end_run()), this is undefined, so I can tell when
# we are in a run or not.
our $basedir=undef;
# Also, the current run number.
our $run_num=-1;


# Max number of sockets to open while programming boxes:
my $max_socket=5;

# The port to connect to on boxes for other control:
my $control_port=23;
# The UDP port for the UDP ping, and message to send (as well as the
# expected response).
my $box_udp_port=11235;
my $box_udp_msg="Go to blackout!";
my $box_udp_rsp="Goin' dark.";
# The TCP and UDP protocols:
my $tcp_protocol=getprotobyname('tcp') || croak("Could not get TCP protocol ($!).\n");
my $udp_protocol=getprotobyname('udp') || croak("Could not get UDP protocol ($!).\n");
# A timeout for waiting for requests.
my $sock_timeout=5;
# The name of the master box -- the one that triggers the rest of the
# boxes (currently, it needs to be an analog box):
my $master_box='analog1';


# A socket to connect to the camera server, and the name of the socket to
# use.  Also, the PID of the camera server, for sending signals.
my $cam_sock=undef;
my $cam_PID=undef;
my $cam_socket='/tmp/camserver';


# This hash contains connections to imgserver programs, for displaying
# images.  The keys are the names of the windows (also included is a
# default name).  The values are command structures for the Fork_System.
my %img_comm=();
my $default_img_name="Experiment Picture";
# Name of program:
my $img_prog='imgserver';
# The task queue for the imgservers.  Allow for many imgservers.
my $img_fork=Fork_System->new_fork();
$img_fork->verbosity(0);
$img_fork->max(10);


# The ping program.
my $ping_prog='/sbin/ping';
# A fork system for the ping program.
my $ping_fork=Fork_System->new_fork();
$ping_fork->verbosity(0);
$ping_fork->max(10);
# The pwd program.
my $pwd_prog='/bin/pwd';


# The fork structure for the user-run commands.
my $user_fork=Fork_System->new_fork();
$user_fork->verbosity(0);
$user_fork->max(5);


# REAP ROUTINE ########################
# If you take over SIGCHLD, you need to call this with the PID and return
# value of every child that may have been spawned by this routine (it
# ignores things that do not belong to it).

sub reap($$){
	my($PID,$ret)=@_;
	# Ignore errors.
	$img_fork->reap($PID,$ret);
	$ping_fork->reap($PID,$ret);
	$user_fork->reap($PID,$ret);
}
# Trap SIGCHLD -- just call my reap() routine.
$SIG{'CHLD'}=sub{
	my $p;while(($p=waitpid(-1,WNOHANG))>0){reap($p,$?);}
};


# KILL ROUTINE ########################
# Very important -- when you hit CTRL-C, be sure to kill all the children.
# Otherwise, things could get annoying.
sub killchildren($){
	# Make sure the forks don't start anything new:
	$img_fork->max(0);$ping_fork->max(0);$user_fork->max(0);
	# Send a kill signal (not a hard kill):
	$img_fork->killall();$ping_fork->killall();$user_fork->killall();
	# Let them die on their own -- don't wait for them, or follow up with a
	# hard kill.
	print "Received signal.  Quitting.\n";exit(@_);
}
$SIG{'INT'}=$SIG{'QUIT'}=$SIG{'TERM'}=\&killchildren;


# INTERNAL ROUTINES ###################


# Attempts to open a socket to a given (host,port) combination (TCP).
# Prints error messages.
# Returns the socket on success, and undef on failure.
sub open_socket($$){
	my($host,$port)=@_;
	my($sock,$iaddr,$paddr,$proto,$timestruct);

	# Some checks.
	if(!$host){carp("open_socket(): bad hostname.\n");return undef;}
	if(!$port){carp("open_socket(): bad port.\n");return undef;}
	if($port=~/\D/){
		if(!($port=getservbyname($port,"tcp"))){
			carp("open_socket(): unable to get port number ($!).\n");
			return undef;
		}
	}

	# Get some important stuff.
	if(!($iaddr=inet_aton($host))){
		carp("open_socket(): unable to create address ($!).\n");
		return undef;
	}
	if(!($paddr=sockaddr_in($port,$iaddr))){
		carp("open_socket(): unable to create socket address ($!).\n");
		return undef;
	}

	# Make the socket and connection.
	if(!socket($sock,,PF_INET,SOCK_STREAM,$tcp_protocol)){
		carp("open_socket(): unable to create socket ($!).\n");
		return undef;
	}
	{
		my($sec,$microsec);
		$sec=POSIX::floor($sock_timeout);$microsec=($sock_timeout-$sec)*1000000;
		$timestruct=pack('LL',$sec,$microsec);
	}
	if(!setsockopt($sock,SOL_SOCKET,SO_SNDTIMEO,$timestruct) ||
	   !setsockopt($sock,SOL_SOCKET,SO_RCVTIMEO,$timestruct)){
		carp("open_socket(): unable to set socket timeouts ($!).\n");
		return undef;
	}
	# The timeouts above don't seem to apply to connect, so do it this way:
	{
		my $flag;
		if(!($flag=fcntl($sock,F_GETFL,0))){
			carp("open_socket(): unable to fetch socket flags ($!).\n");
			return undef;
		}
		if(!fcntl($sock,F_SETFL,$flag | O_NONBLOCK)){
			carp("open_socket(): unable to set socket flags ($!).\n");
			return undef;
		}
		if(!connect($sock,$paddr) && ($! != EINPROGRESS)){
			carp("open_socket(): unable to connect socket ($!).\n");
			return undef;
		} 
		if(!fcntl($sock,F_SETFL,$flag)){
			carp("open_socket(): unable to set socket flags ($!).\n");
			return undef;
		}
	}
	# The connect call was done in non-blocking mode, so now, wait the
	# appropriate amount.
	{
		my $end_time=Time::HiRes::time()+$sock_timeout;
		my($bits,$time);
		while(($time=Time::HiRes::time())<$end_time){
			$bits='';vec($bits,fileno($sock),1)=1;
			select(undef,$bits,undef,$end_time-$time);
			if(vec($bits,fileno($sock),1)){last;}
		}
		if($time>=$end_time){
			carp("open_socket(): timed out waiting for connection.\n");
			return undef;
		}
	}
	# Done.  Make it binary and autoflushing, and then return.
	my $oldFH=select($sock);$|=1;select($oldFH);
	binmode($sock);return $sock;
}


# Sends a string over the given socket.
# Don't buffer anything, so use syswrite.
# Arguments are socket handle and string to send.
# Returns bytes written for success, undef for a write error.
sub write_socket($$){
	my($sock,$strg)=@_;
	my($s,$size);
	if(!$strg || !length($strg)){return 0;}
	# Check for an open connection, and whether $sock seems to be open.
	if(!defined(fileno($sock))){return undef;}
	$size=0;for(;;){
		# Check if the connection closed on us or not (don't check the write
		# flags -- they don't give me much useful information):
		my $bits='';vec($bits,fileno($sock),1)=1;
		select(undef,undef,$bits,0);
		if(vec($bits,fileno($sock),1)){return undef;}
		$s=syswrite($sock,substr($strg,$size));
		if(!$s){last;}
		$size+=$s;($size>=length($strg)) && last;
	}
	($size<length($strg)) && return undef;
	return $size;
}


# Reads data from the connection, unbuffered.
# Returns what was read, which will be an empty string if nothing happened,
# and undef if there was an error.  A likely error is that the connection
# was closed from the other end.  The only errors are if the socket is not
# open or the connection was closed.
sub read_socket($){
	my($sock)=@_;
	my($bits,$ebits,$strg,$end_time)=('','','',Time::HiRes::time()+$sock_timeout);
	# Check for an open connection, and whether $sock seems to be open.
	if(!defined(fileno($sock))){return undef;}
	for(;;){
		my $time=Time::HiRes::time();
		if($time>$end_time){last;}
		$bits='';$ebits='';
		vec($bits,fileno($sock),1)=1;
		vec($ebits,fileno($sock),1)=1;
		select($bits,undef,$ebits,$end_time-$time);
		if(vec($ebits,fileno($sock),1)){return undef;}
		if(vec($bits,fileno($sock),1)){
			sysread($sock,$strg,65536);
			# If select says there was something to read, but there wasn't, that
			# seems to be an indication that the connection was closed.
			if(length($strg)==0){return undef;}
			return $strg;
		}
	}
}


# A simple little routine that sends data out a socket and then waits for
# the given response.  Arguments are the socket handle, string to send, and
# the expected reponse.
# Returns undef for error, 1 for successful response, and 0 for no matching
# response.
sub expect($$$){
	my($sock,$send,$expect)=@_;
	my($response,$r,$end_time)=('',0,Time::HiRes::time()+$sock_timeout);
	# Only send if there is something to send:
	if(defined($send) && length($send)){
		if(!write_socket($sock,$send."\r\n")){return undef;}
	}
	for(;;){
		my $time=Time::HiRes::time();
		# Quit on timeout or error.  Otherwise, add it to what I have and check
		# for a match.
		if($time>$end_time){last;}
		$r=read_socket($sock);if(!defined($r)){return undef;}
		if(!length($r)){last;}$response.=$r;
		if($response=~/$expect/){return 1;}
	}
	return 0;
}


# Basically just a macro to convert an address to a host and port string.
# This is used as a key to the hash of what hosts to ping as well as for
# printing, so it needs to be unique and have everything it needs in it.
# This is used as a hash key, because the actual address structure for a
# response may not be identical to the one for what was sent.
# Returns undef for error, and a string for success.
sub addr($){
	my($paddr)=@_;
	if(!defined($paddr)){return undef;}
	my($port,$iaddr)=sockaddr_in($paddr);
	if(!defined($port) || !defined($iaddr)){return undef;}
	my $host=inet_ntoa($iaddr);if(!defined($host)){return undef;}
	return "$host:$port";
}


# Given a file handle, reads everything it can up until there is no more to
# read -- data is discarded.  Returns true normally, and 0 if the handle
# appears to be closed.
sub dump_handle($){
	my($hand)=@_;
	my($bits,$data);
	for(;;){
		$bits='';vec($bits,fileno($hand),1)=1;
		select($bits,undef,undef,0);
		if(!(vec($bits,fileno($hand),1))){return 1;}
		sysread($hand,$data,65536);
		if(!length($data)){return 0;}
	}
}


# Given a hash ref, this treats each key as a host name, and sends a UDP
# ping to each one.  It then waits for responses until it has waited for a
# specified time, at which point it returns.  The return value is undef for
# error, and a hash ref for success (the keys of the hash are the hosts
# that responded, and the values are false for bad responses and true for
# valid responses).  Complains in the event of an error.
sub udp_ping($){
	my($ref)=@_;
	my($box,$addr,$sock,%rsp,%ret);
	# Checks:
	if(ref($ref) ne 'HASH'){
		carp("udp_ping(): bad argument.\n");
		return undef;
	}
	# Make socket:
	if(!socket($sock,PF_INET,SOCK_DGRAM,$udp_protocol)){
		carp("udp_ping(): error creating socket ($!).\n");
		return undef;
	}
	binmode($sock);
	# Ping:
	foreach $box (keys(%{$ref})){
		my($size,$s);
		# Make an address:
		$addr=inet_aton($box);if(!defined($addr)){
			carp("udp_ping(): unable to make an address for $box.\n");
			return undef;
		}
		$addr=sockaddr_in($box_udp_port,$addr);if(!defined($addr)){
			carp("udp_ping(): unable to make socket address for $box.\n");
			return undef;
		}
		$rsp{addr($addr)}={'box'=>$box};
		# Send the ping:
		for($size=0;;){
			$s=send($sock,substr($box_udp_msg,$size),0,$addr);
			if(!defined($s)){
				carp("udp_ping(): error sending ping to $box ($!).\n");
				return undef;
			}
			$size+=$s;if($size>=length($box_udp_msg)){last;}
		}
	}
	# Wait for responses.
	%ret=();
	my($end_time,$time);$end_time=Time::HiRes::time()+$sock_timeout;
	# As long as there are non-responders ...
	while(scalar(keys(%rsp)) && ($time=Time::HiRes::time())<$end_time){
		my($bits,$data,$done);
		$bits='';vec($bits,fileno($sock),1)=1;
		select($bits,undef,undef,$end_time-$time);
		if(vec($bits,fileno($sock),1)){
			# A response!
			# recv() is a little nicer than sysread() for UDP, because you get
			# the address, too.
			$addr=recv($sock,$data,length($box_udp_rsp),0);
			if(!defined($addr)){
				carp("udp_ping(): error reading socket ($!).\n");
				return undef;
			}else{
				$addr=addr($addr);
				# Add the response to the stuff we've gotten, and check to see if
				# it is valid.
				if(defined($rsp{$addr})){
					if(!defined($rsp{$addr}{'rsp'})){$rsp{$addr}{'rsp'}=$data;}
					else{$rsp{$addr}{'rsp'}.=$data;}
					if($rsp{$addr}{'rsp'} eq $box_udp_rsp){
						# Mark a valid response and remove from our wait list.
						$ret{$rsp{$addr}{'box'}}=1;delete($rsp{$addr});
					}elsif(length($rsp{$addr}{'rsp'}) >= length($box_udp_rsp)){
						# This response isn't going to be right -- mark it as wrong,
						# and remove it from our wait list.
						$ret{$rsp{$addr}{'box'}}=0;delete($rsp{$addr});
					}
				}
			}
		# End of if(received response){...}
		}
	# End of loop waiting for reponses.
	}
	# Mark the rest of the bad-responders.
	foreach $addr (keys(%rsp)){
		if(defined($rsp{$addr}{'rsp'})){
			$ret{$rsp{$addr}{'box'}}=0;
		}
	}
	return \%ret;
}


# Given a list of actions (described below -- these had better be sorted by
# time), this forms a hash.  The keys are all the digital boxes referenced
# (by name), and the values are hash references, of which the important
# keys are 'data', which contain the string with which the box should be
# programmed.  A more complete description of the hash is in the comments
# preceding the program_boxes() routine (the box_data data type).
# Returns the hash reference if successful and undef for error.  
# You can go ahead and give all the actions here -- it only pays attention
# to the Digital_Action s.
# Warns about various likely-to-be-problems.
# When conflicting actions occur in the list at the same time, the last one
# is used (as though it comes just a tiny bit after the earlier one in the
# list).
# This is a 'prog_data' routine (see description of %Action near the top of
# this file).  It also needs to set a 'UDP' element in each box's hash,
# since all boxes that are to be programmed should be pinged.
sub digital_box_data(@){
	my @action_list=@_;
	my($action,$box,$delay,$bit,$b);
	my %digital_box_data=();
	foreach $action (@action_list){
		# Ignore things not meant for this routine.
		if(ref($action) ne "Digital_Action"){next;}
		# Have the time in the proper units for this type of box.
		# Although the two boards are 180 degrees out of phase, each gets
		# updated every 4 microseconds (they run off a 10MHz clock, divided
		# down by 10, which in turn runs a flip-flop, and then it alternates
		# between the two boards, for a frequency of 10MHz/40 = 250kHz).
		# FIXME -- handle different times for different boards inside the box?
		my $dig_time=POSIX::floor($action->{'time'}/0.004+0.5);
		$box=$action->{'box'};
		if(!defined($digital_box_data{$box})){
			$digital_box_data{$box}={
				'action_type'=>'Digital_Action',
				'data'=>'',
				'time'=>$dig_time,
				'state'=>0x00000000,
				'new_hi'=>0x00000000,
				'new_lo'=>0x00000000,
				# We want this box to be UDP pinged:
				'UDP'=>1,
			};
		}
		# If the current time is after the current update for the box, add
		# the data and update the time (but don't add the data for negative
		# time -- that's just setting the initial condition).
		# FIXME -- update only if state has changed?  or not, so you can have
		#   delays longer than max time?
		if($dig_time>$digital_box_data{$box}{'time'}){
			# A delay of 0 means the update will be on the next cycle, so we
			# are off by one (hence the "-1").
			# Also, make sure that, as far as programming goes, time actually
			# starts at 0, not whatever negative event last set the state.
			# We won't be making a real change until $dig_time is past zero,
			# not just at 0.
			if($dig_time>0 && $digital_box_data{$box}{'time'}<0){
				$digital_box_data{$box}{'time'}=0;
			}
			$delay=$dig_time-$digital_box_data{$box}{'time'}-1;
			if($delay>4000000000){
				carp("digital_box_data(): obscenely large delay ($delay) at time ",$action->{'time'}," for $box.\n");
				return undef;
			}
			$digital_box_data{$box}{'state'}|=( $digital_box_data{$box}{'new_hi'});
			$digital_box_data{$box}{'state'}&=(~$digital_box_data{$box}{'new_lo'});
			# Don't actually write a change until the first time you should
			# change it -- AFTER time 0.
			if($dig_time>0){
				$digital_box_data{$box}{'data'}.=
					pack("VV",($delay,$digital_box_data{$box}{'state'}));
			}
			$digital_box_data{$box}{'new_hi'}=0x00000000;
			$digital_box_data{$box}{'new_lo'}=0x00000000;
			$digital_box_data{$box}{'time'}=$dig_time;
		}
		# Handle the event.
		foreach $bit (@{$action->{'hi'}}){
			$b=(1<<$bit);
			if($digital_box_data{$box}{'new_lo'} & $b){
				carp("digital_box_data(): WARNING: attempting to set bit $bit both hi and lo at time ",$action->{'time'}," on $box.\n");
				$digital_box_data{$box}{'new_lo'}&=(~$b);
			}
			if($digital_box_data{$box}{'new_hi'} & $b){
				carp("digital_box_data(): WARNING: trying multiple times to set bit $bit hi at time ",$action->{'time'}," on $box.\n");
			}else{$digital_box_data{$box}{'new_hi'} |= $b;}
		}
		foreach $bit (@{$action->{'lo'}}){
			$b=(1<<$bit);
			if($digital_box_data{$box}{'new_hi'} & $b){
				carp("digital_box_data(): WARNING: attempting to set bit $bit both hi and lo at time ",$action->{'time'}," on $box.\n");
				$digital_box_data{$box}{'new_hi'}&=(~$b);
			}
			if($digital_box_data{$box}{'new_lo'} & $b){
				carp("digital_box_data(): WARNING: trying multiple times to set bit $bit lo at time ",$action->{'time'}," on $box.\n");
			}else{$digital_box_data{$box}{'new_lo'} |= $b;}
		}
	# End of digital box programming loop for a given action.
	}
	# Write whatever else needs to be written, to set the final state, and
	# clean up extra stuff that will not be used again.
	foreach $box (keys(%digital_box_data)){
		if($digital_box_data{$box}{'new_hi'} ||
		   $digital_box_data{$box}{'new_lo'}){
			$digital_box_data{$box}{'state'}|=( $digital_box_data{$box}{'new_hi'});
			$digital_box_data{$box}{'state'}&=(~$digital_box_data{$box}{'new_lo'});
			$digital_box_data{$box}{'data'}.=
				pack("VV",(0,$digital_box_data{$box}{'state'}));
			$digital_box_data{$box}{'new_hi'}=0x00000000;
			$digital_box_data{$box}{'new_lo'}=0x00000000;
		}
		delete($digital_box_data{$box}{'time'});
		delete($digital_box_data{$box}{'state'});
		delete($digital_box_data{$box}{'new_hi'});
		delete($digital_box_data{$box}{'new_lo'});
	}
	return \%digital_box_data;
}


# Given a list of actions (described below -- these had better be sorted by
# time), this forms a hash.  The keys are all the digital boxes referenced
# (by name), and the values are hash references, of which the important
# keys are 'data', 'idata', and 'init', which contain the string with which
# the box should be programmed and data for the initial setup (and
# something to save to a file to record what the initial setup was).  A
# more complete description of the hash is in the comments preceding the
# program_boxes() routine (the box_data data type).  Returns the hash
# reference if successful, and undef for error.  
# You can go ahead and give all the actions here -- it only pays attention
# to the Analog_Action s.
# Warns about various likely-to-be-problems.
# When conflicting actions occur in the list at the same time, the last one
# is used (as though it comes just a tiny bit after the earlier one in the
# list).
# This is a 'prog_data' routine (see description of %Action near the top of
# this file).  It also needs to set a 'UDP' element in each box's hash,
# since all boxes that are to be programmed should be pinged.
sub analog_box_data(@){
	my @action_list=@_;
	my($action,$box);
	my %analog_box_data=();
	# FIXME -- assuming a fixed update time for the boards:
	my $update_time=0.004;
	# FIXME -- assuming a fixed update time for the boards, this is the
	# number of time steps I consider to be "too long":
	my $long_time=4;
	# This first loop is to sort the actions based on which box they affect,
	# and within each box, based on the direct sets, which have the highest
	# priorities, and the interpolated ramps, which will be used to fill in
	# the gaps afterwards.
	foreach $action (@action_list){
		# Ignore things not meant for this routine.
		if(ref($action) ne "Analog_Action"){next;}
		# Have the time in the proper units for this type of box.
		# You can only update one output at a time, and the update rate can be
		# set on each board.
		# FIXME -- assuming a fixed update time for the boards.
		my $ana_time=POSIX::floor($action->{'time'}/$update_time+0.5);
		$box=$action->{'box'};
		if(!defined($analog_box_data{$box})){
			$analog_box_data{$box}={
				'action_type'=>'Analog_Action',
				# Two lists -- one for the direct set actions, and one for the
				# interpolations.
				'set' =>[],
				'ramp'=>[],
				# We want this box to be pinged:
				'UDP'=>1,
			};
		}
		# Add the action to the right list, expanding ramps.
		if($action->{'type'} eq 'set'){
			my %h=(
				'time'=>$ana_time,
				'orig'=>$ana_time,
				'ord'=>defined($action->{'ord'})?$action->{'ord'}:0,
				'out'=>$action->{'out'},
				'val'=>$action->{'val'},
			);
			push(@{$analog_box_data{$box}{'set'}},\%h);
		}elsif($action->{'type'} eq 'ramp'){
			# First, find when the last change was.  If there wasn't a last
			# change, make up an initial value.
			my($init_t,$init_v)=(0,0);
			# I'm going to guess this is more efficient than a for loop:
			foreach my $a (reverse(@{$analog_box_data{$box}{'set'}})){
				if($a->{'out'}==$action->{'out'}){
					$init_t=$a->{'time'};
					$init_v=$a->{'val'};
					last;
				}
			}
			# The loop:  don't bother with a bunch of negative times, and also
			# make sure that the last one is a set point, so that at least the
			# end of the ramp gets some priority.  Also, don't need to reset the
			# initial point.  Get one negative time, so if a ramp starts
			# somewhere negative, I can initialize the channel based on where the
			# ramp should be at time=0 (or, -1).
			my $last_v=POSIX::floor(($init_v+10)*0x10000/20+0.5);
			my $goal_v=POSIX::floor(($action->{'val'}+10)*0x10000/20+0.5);
			if($last_v>0xFFFF){$last_v=0xFFFF;}if($goal_v>0xFFFF){$goal_v=0xFFFF;}
			# It turns out that very slow ramps lug perl down quite a bit -- a
			# 1 second ramp has 250,000 time steps, and 10 seconds breaks a
			# million.  Perl takes a long time to go through those loops, so
			# detect the cases where time changes more rapidly than value, and
			# use a different loop.
			if(abs($goal_v-$last_v)>=$ana_time-$init_t){
				for(my $t=($init_t<-1)?-1:$init_t+1;$t<$ana_time;++$t){
					my %h=(
						'time'=>$t,
						'orig'=>$t,
						# For repeatability, try to make ord unique based on the ord of
						# the event that caused this.
						'ord'=>defined($action->{'ord'})?(-$action->{'ord'}):-1,
						'out'=>$action->{'out'},
						# Note that the for loop guarantees that
						# $action->{'time'}!=$init_t:
						'val'=>$init_v+($action->{'val'}-$init_v)*($t-$init_t)/($ana_time-$init_t),
					);
					# Only add it if there is a change in the output.
					my $val=POSIX::floor(($h{'val'}+10)*0x10000/20+0.5);
					if($val>0xFFFF){$val=0xFFFF;}
					if($val!=$last_v){
						push(@{$analog_box_data{$box}{'ramp'}},\%h);
						$last_v=$val;
					}
				}
			}else{
				# Here is the version where time changes more often than values:
				# Note that it first and last points are placed right at the times
				# they are supposed to be placed, but that other points are placed
				# as close to a line as I can get.  The deviation from linearity at
				# the endpoints should be small.
				if($goal_v<$last_v){
					# The for loop bounds, and the fact that $last_v and $goal_v
					# should be integers should deal with possible divide-by-0
					# problems.
					for(my $v=$last_v-1;$v>$goal_v;--$v){
						my $t=POSIX::floor($init_t+($v-$last_v)*($ana_time-$init_t)/($goal_v-$last_v)+0.5);
						# More efficient to jump to where $t>=-1, I suppose, but this is
						# so easy:
						if($t<-1){next;}
						my %h=(
							'time'=>$t,
							'orig'=>$t,
							# For repeatability, try to make ord unique based on the ord of
							# the event that caused this.
							'ord'=>defined($action->{'ord'})?(-$action->{'ord'}):-1,
							'out'=>$action->{'out'},
							# I already have the scaled version, so store that somewhere
							# special:
							'valint'=>$v,
						);
						# The output is guaranteed to be different from the previous
						# value:
						push(@{$analog_box_data{$box}{'ramp'}},\%h);
					}
				}else{
					for(my $v=$last_v+1;$v<$goal_v;++$v){
						my $t=POSIX::floor($init_t+($v-$last_v)*($ana_time-$init_t)/($goal_v-$last_v)+0.5);
						# More efficient to jump to where $t>=-1, I suppose, but this is
						# so easy:
						if($t<-1){next;}
						my %h=(
							'time'=>$t,
							'orig'=>$t,
							# For repeatability, try to make ord unique based on the ord of
							# the event that caused this.
							'ord'=>defined($action->{'ord'})?(-$action->{'ord'}):-1,
							'out'=>$action->{'out'},
							# I already have the scaled version, so store that somewhere
							# special:
							'valint'=>$v,
						);
						# The output is guaranteed to be different from the previous
						# value:
						push(@{$analog_box_data{$box}{'ramp'}},\%h);
					}
				}
			}
			# Now, the final point.  Give it a higher 'ord', because I don't want
			# it to get overwritten:
			my %h=(
				'time'=>$ana_time,
				'orig'=>$ana_time,
				'ord'=>defined($action->{'ord'})?$action->{'ord'}:0,
				'out'=>$action->{'out'},
				'val'=>$action->{'val'},
			);
			push(@{$analog_box_data{$box}{'set'}},\%h);
		}else{
			carp("analog_box_data(): unknown type ",$action->{'type'}," of Analog_Action for $box.\n");
			return undef;
		}
	# End of first loop through actions.
	}

	# Everything else may now be done on a per box basis:
	foreach $box (keys(%analog_box_data)){
		# Key is a channel; value refers to the last element of the
		# set array that set that channel:
		my %last_chan=();
		my $time=0;
		# Second pass -- run through every element in the set array, removing
		# conflicts by shifting them later in time.  When there's a time conflict
		# in the same channel, issue a warning and overwrite the previous one
		# (they are supposed to be sorted so that later ones take precedence over
		# earlier ones).
		foreach my $a (@{$analog_box_data{$box}{'set'}}){
			# Initializations are handled differently:
			if($a->{'time'}<0){
				# We'll add it to the initialization list, which we'll create if it
				# isn't already made.
				if(!defined($analog_box_data{$box}{'idata'})){
					$analog_box_data{$box}{'idata'}={};
				}
				# Depend on things being sorted in order of precedence -- always
				# overwrite initializations, but warn about it:
				if(defined($analog_box_data{$box}{'idata'}{$a->{'out'}})){
					carp("analog_box_data(): WARNING: multiple initializations of channel ",$a->{'out'}," on $box.\n");
				}
				$analog_box_data{$box}{'idata'}{$a->{'out'}}=$a->{'val'};
				# Don't need this element anymore:
				undef($a);next;
			}
			# Only overwrite a previous channel if both of these really are
			# supposed to go in the same time slot:
			if(defined($last_chan{$a->{'out'}}) &&
			   $last_chan{$a->{'out'}}{'orig'}==$a->{'time'}){
				carp("analog_box_data(): WARNING: multiple sets occuring around time ",$a->{'time'}*$update_time," for channel ",$a->{'out'}," on $box.\n");
				$last_chan{$a->{'out'}}{'ord'}=$a->{'ord'};
				$last_chan{$a->{'out'}}{'val'}=$a->{'val'};
				# $a has now been effectively moved back -- attempt to remove it
				# from the array.
				undef($a);next;
			}
			if($a->{'time'}>=$time){$time=$a->{'time'}+1;}
			else{
				if($time-$a->{'time'}>$long_time){
					carp("analog_box_data(): WARNING: setting output ",$a->{'out'}," to ",$a->{'val'}," on $box is delayed from about ",$a->{'time'}*$update_time," to ",$time*$update_time,".\n");
				}
				$a->{'time'}=$time;++$time;
			}
			# Store in %last_chan:
			$last_chan{$a->{'out'}}=$a;
		}
		# Quickly scan the ramps for negative times -- they override
		# initialization with a warning.
		foreach $a (@{$analog_box_data{$box}{'ramp'}}){
			if($a->{'time'}<0){
				# Warn:
				carp("analog_box_data(): WARNING: you have a ramp at negative time on channel ",$a->{'out'}," on $box.  This overrides other initializations for this channel -- did you intend that?\n");
				# We'll add it to the initialization list, which we'll create if it
				# isn't already made.
				if(!defined($analog_box_data{$box}{'idata'})){
					$analog_box_data{$box}{'idata'}={};
				}
				# Some ramps may have 'valint' rather than 'val'.  That needs to be
				# converted back to a 'val' (and, during initialization, will be
				# converted to a voltage on the Ethernut running the analog box ...
				# lots of converting!).
				if(defined($a->{'val'})){
					$analog_box_data{$box}{'idata'}{$a->{'out'}}=$a->{'val'};
				}elsif(defined($a->{'valint'})){
					$analog_box_data{$box}{'idata'}{$a->{'out'}}=$a->{'valint'}*20/0x10000-10;
				}
			}
		}
		# Third pass:  add the ramp array elements to the set array.  This
		# differs from before in that there is no set order for the stuff.
		# Look for optimum places to put what you can.  Look both ahead and
		# before the given time.  If you run into a set of the same channel,
		# consider replacing it.
		# I'm going to want to traverse the set array quickly and easily,
		# inserting elements willy-nilly, so change it to a linked list.
		# IMPORTANT:  This is a bi-directional linked list.  That means every
		# element points to the next and previous element.  Thus, we have all
		# these little elements that are anonymous hashes, and they contain
		# references to each other.  Even when the subroutine ends, and the
		# reference to all of these go out of scope, because they all have
		# references to each other (the 'next' and 'prev' guys), perl cannot
		# deallocate these individual hashes.  Thus, as you keep calling this,
		# these guys build up.  If you repeatedly call this with hefty ramps,
		# you can fill up the memory on your computer.
		# THE ABOVE RESULTS IN A MEMORY LEAK.
		# The fix: manually remove the references whenever the subroutine
		# returns.
		# This is mentioned in the Two-Phased Garbage Collection section of the
		# perlobj man page.
		my($linked_set,$last_link);
		foreach my $s (@{$analog_box_data{$box}{'set'}}){
			# Skip the undefined ones.
			if(!defined($s)){next;}
			if(!defined($linked_set)){
				$linked_set={'ref'=>$s};$last_link=$linked_set;
			}else{
				my %temp=('ref'=>$s,'prev'=>$last_link);
				$last_link->{'next'}=\%temp;$last_link=\%temp;
			}
		}
		# Don't need this anymore:
		delete($analog_box_data{$box}{'set'});
		# The ramp array needs sorting.
		my @ramp=sort {
			my $ret=($a->{'time'} <=> $b->{'time'});if($ret){return $ret;}
			# Do the ord in reverse, so that smaller priorities get stuck AFTER
			# lower priorities.
			return($b->{'ord'} <=> $a->{'ord'});
		} @{$analog_box_data{$box}{'ramp'}};
		# Don't need this anymore:
		delete($analog_box_data{$box}{'ramp'});
		# NOTE:  Allow for the case that $linked_set is not defined.
		my $s=$linked_set;
		while(scalar(@ramp)){
			# Collect all elements with the same time.
			my @r=(shift(@ramp));
			while(scalar(@ramp)){
				if($ramp[0]{'time'}!=$r[0]{'time'}){last;}
				push(@r,shift(@ramp));
			}
			# Skip the negative-time ramps -- they were for initialization:
			if($r[0]{'time'}<0){next;}
			# Find where we are in the set list.
			if(defined($s)){
				while($s->{'ref'}{'time'}<$r[0]{'time'}){$s=$s->{'next'};}
			}
			# @r now has all the elements competing for the same time.
			# Sort them based on how long ago each channel was set, or how long
			# For speed, only bother with this if there is a conflict (more than
			# one trying to get this time slot):
			if(scalar(@r)>1){
				# until it will next be set, whichever is shorter.
				my %set_age=();
				# Find the set age for each channel represented in @r -- not the
				# fastest way to do it, but not too complicated.
				foreach my $o (@r){
					my $c;
					# In the ramp array, there shouldn't ever be multiple sets for the
					# same channel at a given time.
					if(defined($set_age{$o->{'out'}})){
						carp("analog_box_data(): Internal bug: multiple sets in ramp array for single channel while computing ramps for $box.\n");
						# Avoid evil memory leaks by removing self-reference loops in
						# the linked list (see above, where it says MEMORY LEAK):
						my $s=$linked_set;while(defined($s)){
							my $n=$s->{'next'};
							delete($s->{'prev'});delete($s->{'next'});
							$s=$n;
						}
						return undef;
					}
					if(defined($s)){
						# Look forwards in linked_set:
						$c=$s;while($c->{'ref'}{'out'}!=$o->{'out'}){
							if(!defined($c->{'next'})){last;}
							$c=$c->{'next'}
						}
						if($c->{'ref'}{'out'}==$o->{'out'}){
							if(!defined($set_age{$o->{'out'}}) ||
							   $c->{'ref'}{'time'}-$o->{'time'}<$set_age{$o->{'out'}}){
								$set_age{$o->{'out'}}=$c->{'ref'}{'time'}-$o->{'time'};
							}
						}
						# Look backwards in linked_set:
						$c=$s;while($c->{'ref'}{'out'}!=$o->{'out'}){
							if(!defined($c->{'prev'})){last;}
							$c=$c->{'prev'}
						}
						if($c->{'ref'}{'out'}==$o->{'out'}){
							if(!defined($set_age{$o->{'out'}}) ||
							   $o->{'time'}-$c->{'ref'}{'time'}<$set_age{$o->{'out'}}){
								$set_age{$o->{'out'}}=$o->{'time'}-$c->{'ref'}{'time'};
							}
						}
					}
					# Look forwards in ramp:
					foreach $c (@ramp){
						if($c->{'out'}==$o->{'out'} &&
						   $c->{'time'}-$o->{'time'}<$set_age{$o->{'out'}}){
							$set_age{$o->{'out'}}=$c->{'time'}-$o->{'time'};
							last;
						}
					}
				}
				# Now that we have the set age for each one, sort them based on it:
				@r=sort {
					if(!defined($set_age{$a->{'out'}})){
						if(!defined($set_age{$b->{'out'}})){return 0;}
						return 1;
					}
					if(!defined($set_age{$b->{'out'}})){return -1;}
					my $ret=($set_age{$b->{'out'}} <=> $set_age{$a->{'out'}});
					if($ret){return $ret;}
					return($b->{'ord'} <=> $a->{'ord'});
				} @r;
			# End of more than one ramp element competing for this time:
			}
			# Now add each element of @r to my linked_set -- look both
			# directions:
			my($p,$n);while(scalar(@r)){
				my $r=pop(@r);
				my $c;
				# If $s is undefined, then the set array is currently empty, so I
				# can just add this wherever I want.
				if(!defined($s)){
					$s={'ref'=>$r};
					$linked_set=$s;next;
				}
				# Reset these, because I need to search for a similar output, in
				# which case I can delete this.
				$p=$n=$s;
				# Go forwards to find a spot, and set $n to the one to follow (n
				# for next) ($n should have a time>=what I want, and previous
				# increments should only have skipped things with no room to spare):
				while(defined($n->{'next'})){
					$c=$n->{'next'};
					if($c->{'ref'}{'time'}-$n->{'ref'}{'time'}>1){last;}
					if($n->{'ref'}{'out'}==$r->{'out'}){last;}
					$n=$c;
				}
				# Go backwards to find another spot, and set $p to the one to
				# precede (p for precede):
				while(defined($p->{'prev'})){
					$c=$p->{'prev'};
					if($p->{'ref'}{'time'}-$c->{'ref'}{'time'}>1){last;}
					if($p->{'ref'}{'out'}==$r->{'out'}){last;}
					$p=$c;
				}
				# Pick where to put it.  $where=1 means try to replace $n, 2 means
				# after $n, -1 means replace $p, and -2 means before $p.  0 is
				# right here.
				my $where=undef;
				# If it can go at the current time, then $p must have a time
				# greater than the current time, and we can safely insert before
				# $p:
				if($p->{'ref'}{'time'}>$r->{'time'}){$where=0;}
				else{
					# If $p points to a previous change of this channel ...
					if($p->{'ref'}{'out'}==$r->{'out'}){
						# ... and $n does not, choose after $n:
						if($n->{'ref'}{'out'}!=$r->{'out'}){$where=2;}
						# Otherwise, replace the closest one:
						elsif($r->{'time'}-$p->{'ref'}{'time'}>$n->{'ref'}{'time'}-$r->{'time'}){
							$where=1;
						}else{$where=-1;}
					}else{
						# ... if $n points to a change of this channel, go before $p
						# ...
						if($n->{'ref'}{'out'}==$r->{'out'}){
							# ... or replace $n if before $p is not an option.
							if($p->{'ref'}{'time'}<=0){$where=1;}
							else{$where=-2;}
						}
						# Otherwise, go before or after the closest one (or just after
						# if before $p is not an option):
						elsif($r->{'time'}-$p->{'ref'}{'time'}>$n->{'ref'}{'time'}-$r->{'time'}){
							$where=2;
						}elsif($p->{'ref'}{'time'}<=0){$where=2;}
						else{$where=-2;}
					}
				}
				# Shouldn't happen:
				if(!defined($where)){
					carp("analog_box_output(): Internal bug: \$where is undefined on $box, out ",$r->{'out'},".\n");
					# Avoid evil memory leaks by removing self-reference loops in
					# the linked list (see above, where it says MEMORY LEAK):
					my $s=$linked_set;while(defined($s)){
						my $n=$s->{'next'};
						delete($s->{'prev'});delete($s->{'next'});
						$s=$n;
					}
					return undef;
				}
				# Now, put it away:
				if($where==0 || $where==-2){
					# If we really are going before $p, set the time.
					if($where==-2){$r->{'time'}=$p->{'ref'}{'time'}-1;}
					# If we really are putting it right here, we may need to put it
					# forwards a notch or two.  There's definitely a gap before $p,
					# but the beginning of the gap might begin right here at
					# $r->{'time'} (if we've put one of the earlier competitors for
					# this time slot there -- that wouldn't move $s, so $p->{'prev'}
					# would have a time equal to the current time).
					if($where==0 && defined($p->{'prev'}) &&
					   $p->{'prev'}{'ref'}{'time'}>=$r->{'time'}){
						$r->{'time'}=$p->{'prev'}{'ref'}{'time'}+1;
					}
					my %temp=('ref'=>$r,'next'=>$p);
					# In theory, we are adding things that occur AFTER a set on the
					# channel, so we should never be adding to the beginning of the
					# list, since we won't go back past a channel set.  However,
					# there may be some weird thing with non-initialized channels, so
					# allow for it.
					if(defined($p->{'prev'})){
						$temp{'prev'}=$p->{'prev'};$p->{'prev'}{'next'}=\%temp;
					}else{$linked_set=\%temp;}
					$p->{'prev'}=\%temp;
				}elsif($where==2){
					$r->{'time'}=$n->{'ref'}{'time'}+1;
					my %temp=('ref'=>$r,'prev'=>$n);
					if(defined($n->{'next'})){
						$temp{'next'}=$n->{'next'};
						$n->{'next'}{'prev'}=\%temp;
					}
					$n->{'next'}=\%temp;
				}elsif($where==1){
					# Only replace $n if it is NOT a set, and ours is closer to its
					# original time than $n.  Shouldn't be here if the 'out' values
					# do not match.
					if($n->{'ref'}{'ord'}<0 &&
					   abs($n->{'ref'}{'time'}-$n->{'ref'}{'orig'}) >
					   abs($n->{'ref'}{'time'}-$r->{'time'})){
						$n->{'ref'}{'orig'}=$r->{'orig'};
						# Ramp elements may have a valint instead of a val:
						if(defined($r->{'valint'})){$n->{'ref'}{'valint'}=$r->{'valint'};}
						else{$n->{'ref'}{'val'}=$r->{'val'};}
					}
				}elsif($where==-1){
					# Only replace $p if it is NOT a set, and ours is closer to its
					# original time than $p.  Shouldn't be here if the 'out' values
					# do not match.
					if($p->{'ref'}{'ord'}<0 &&
					   abs($p->{'ref'}{'time'}-$p->{'ref'}{'orig'}) >
					   abs($p->{'ref'}{'time'}-$r->{'time'})){
						$p->{'ref'}{'orig'}=$r->{'orig'};
						# Ramp elements may have a valint instead of a val:
						if(defined($r->{'valint'})){$p->{'ref'}{'valint'}=$r->{'valint'};}
						else{$p->{'ref'}{'val'}=$r->{'val'};}
					}
				}else{
					carp("analog_box_output(): Internal bug: \$where is a bad value on $box, out ",$r->{'out'},".\n");
					# Avoid evil memory leaks by removing self-reference loops in
					# the linked list (see above, where it says MEMORY LEAK):
					my $s=$linked_set;while(defined($s)){
						my $n=$s->{'next'};
						delete($s->{'prev'});delete($s->{'next'});
						$s=$n;
					}
					return undef;
				}
			# End of adding @r to $linked_set.
			}
		# End of adding @ramp to $linked_set.
		}

		# Write the initialization string for this box.  Keep it pretty simple:
		# one line per channel that is initialized ("\n" terminated), with
		# channel name (or number), tab, value.
		if(defined($analog_box_data{$box}{'idata'})){
			$analog_box_data{$box}{'init'}='';
			foreach my $a (keys(%{$analog_box_data{$box}{'idata'}})){
				if(defined($analog_box_data{$box}{'idata'}{$a})){
					$analog_box_data{$box}{'init'}.="$a\t$analog_box_data{$box}{'idata'}{$a}\n";
				}
			}
		}
		# Write the program for this box:
		$analog_box_data{$box}{'data'}='';
		# Trivial case -- no program, move on to next box.
		# For now, actually have a null program to send to the box, and still
		# UDP ping it -- why list it if you didn't use it?  Maybe you want the
		# ping.
		if(!defined($linked_set)){next;}
		$s=$linked_set;
		# Write the code for the first bit, if necessary:
		if($s->{'ref'}{'time'}>0){
			# Write to an invalid address:
			$analog_box_data{$box}{'data'}.=
				pack("VCv",$s->{'ref'}{'time'}-1,0xFF,0);
		}
		for(;;){
			# The -1 is because with a delay of 0, the next change actually
			# occurs one cycle later.
			my($n,$delay);
			$n=$s->{'next'};
			if(!defined($n)){$delay=0;}
			else{$delay=$n->{'ref'}{'time'}-$s->{'ref'}{'time'}-1;}
			if($delay<0){
				carp("analog_box_data(): Internal bug: misordered times after ordering.\n");
				# Avoid evil memory leaks by removing self-reference loops in
				# the linked list (see above, where it says MEMORY LEAK):
				my $s=$linked_set;while(defined($s)){
					my $n=$s->{'next'};
					delete($s->{'prev'});delete($s->{'next'});
					$s=$n;
				}
				return undef;
			}
			if($delay>4000000000){
				carp("analog_box_data(): obscenely large delay ($delay) around time ",$s->{'ref'}{'time'}*$update_time,"for $box.\n");
				# Avoid evil memory leaks by removing self-reference loops in
				# the linked list (see above, where it says MEMORY LEAK):
				my $s=$linked_set;while(defined($s)){
					my $n=$s->{'next'};
					delete($s->{'prev'});delete($s->{'next'});
					$s=$n;
				}
				return undef;
			}
			# valint is the already processed version of val.  valint takes
			# priority if it exists.
			my $val;
			if(defined($s->{'ref'}{'valint'})){$val=$s->{'ref'}{'valint'};}
			else{
				$val=POSIX::floor(($s->{'ref'}{'val'}+10)*0x10000/20+0.5);
				if($val>0xFFFF){$val=0xFFFF;}
			}
			if($val<0){
				carp("analog_box_data(): out of range voltage ",$s->{'ref'}{'val'}," on output ",$s->{'ref'}{'out'}," around time ",$s->{'ref'}{'time'}*$update_time,"for $box.\n");
				$val=0;
			}elsif($val>0xFFFF){
				carp("analog_box_data(): out of range voltage ",$s->{'ref'}{'val'}," on output ",$s->{'ref'}{'out'}," around time ",$s->{'ref'}{'time'}*$update_time,"for $box.\n");
				$val=0xFFFF;
			}
			$analog_box_data{$box}{'data'}.=pack("VCv",$delay,$s->{'ref'}{'out'},$val);

			if(!defined($s->{'next'})){last;}
			$s=$s->{'next'};
		}
		# Avoid evil memory leaks by removing self-reference loops in
		# the linked list (see above, where it says MEMORY LEAK):
		while(defined($linked_set)){
			my $n=$linked_set->{'next'};
			delete($linked_set->{'prev'});delete($linked_set->{'next'});
			$linked_set=$n;
		}
	# End of loop for each box.
	}
	return \%analog_box_data;
}


# Given a box_data structure (one particular box, described before the
# program_boxes() routine), this decides if the box should be initialized.
# If not, returns 0.  If so, it initializes the box, returns true for
# success, and undef for error.
# The second argument is the box name, for more useful error messages (and
# so we know which box to which to connect).
sub analog_box_init($$){
	my($box_data,$box)=@_;
	# Check box type:
	if($box_data->{'action_type'} ne 'Analog_Action'){
		carp("analog_box_init(): ERROR: called for wrong box type ".$box_data->{'action_type'}.".\n");
		return undef;
	}
	# Is there any init to be done?
	if(!defined($box_data->{'idata'}) || !scalar(keys(%{$box_data->{'idata'}}))){
		return 0;
	}
	# Initialize:
	# Connect to box:
	my $sock=open_socket($box,$control_port);
	if(!defined($sock)){
		carp("analog_box_init(): unable to open socket to $box.\n");
		return undef;
	}
	# Get to the right menu:
	if(!expect($sock,"","200") ||  # Get the welcome message.
	   !expect($sock,"o","400")){
		carp("analog_box_init(): did not get main menu for box $box.\n");
		close($sock);return undef;
	}
	# Set channels:
	foreach my $a (keys(%{$box_data->{'idata'}})){
		# One trick is that the "telnet thread" uses a different numbering
		# scheme, so translate:
		my $b=16-$a;
		if(!expect($sock,"V","412") ||
		   !expect($sock,$b,"Voltage") ||
		   !expect($sock,$box_data->{'idata'}{$a},"set to")){
			carp("analog_box_init(): error initializing channel ".$a." to ".$box_data->{'idata'}{$a}." on $box.\n");
			close($sock);return undef;
		}
	}
	# Close connection:
	if(!expect($sock,"q","200")){
		carp("analog_box_init(): error getting bach to main menu with $box.\n");
		close($sock);return undef;
	}
	# Don't bother checking for more errors:
	write_socket($sock,"q\r\n");close($sock);
	return 1;
}


# Given a box_data structure (one particular box, described before the
# program_boxes() routine), this interprets the data that was sent, and
# responds to any requests.  The second argument is the most recently read
# stuff from the box.  The third argument is the box name, for more useful
# error messages (this prints error messages).  Returns undef for error,
# negative if box said it was full, 0 if things are still going, and
# postive if we have finished sending data and the connection may be
# closed.  program_boxes() calls this on opening the socket with the second
# argument undefined, in case we need to initiate communications.  We
# don't, so we just return here for that.  In that case, the return value
# is ignored (unless it is undef, which means an error occurred).
# This routine needs to move stuff from the 'data' key to the 'sent' key as
# it is sent.
sub respond_to_box($$$){
	my($box_data,$data,$box)=@_;
	# Return if this was the call that happens when the socket just opens.
	if(!defined($data)){return 0;}
	# Append data to req:
	if(!defined($box_data->{'req'})){$box_data->{'req'}=$data;}
	else{$box_data->{'req'}.=$data;}
	# Now, handle based on box type:
	if($box_data->{'action_type'} eq 'Digital_Action' ||
	   $box_data->{'action_type'} eq 'Analog_Action'){
		# Only respond to a completed request (ends with '\n').
		if($box_data->{'req'}=~/\n$/){
			my $total=0;
			my $req;
			foreach $req (split(/\n+/,$box_data->{'req'})){
				if(!($req=~/^\d+$/)){
					carp("respond_to_box(): bad request from $box.\n");
					next;
				}
				# A request of 0 means it is full.
				if($req==0){return -1;}
				$total+=$req;
			}
			if($total>length($box_data->{'data'})){
				$total=length($box_data->{'data'});
			}
			# $total now has the number of bytes we should send, so send it.
			# For writing to a network, writing should be limited by network
			# speed and stuff, so there's no reason not to let this block.
			# Make a quick check that the socket seems open.
			if(!defined(fileno($box_data->{'sock'}))){
				carp("respond_to_box(): socket for $box does not appear to be open.\n");
				return undef;
			}
			# Ethernuts take data much more smoothly if they get a tiny delay
			# between when they send something, and when you respond.  The easy
			# solution is to insert a tiny delay before sending.  The delay is
			# small enough that it won't noticeable slow down data transfer to
			# other things:
			Time::HiRes::sleep(0.001);
			while(my $s=syswrite($box_data->{'sock'},substr($box_data->{'data'},0,$total))){
				if(!defined($s)){
					carp("respond_to_box(): problem sending data to $box ($!).\n");
					return undef;
				}
				# Move the written portion from data to sent.
				$box_data->{'sent'}.=substr($box_data->{'data'},0,$s);
				$box_data->{'data'}=substr($box_data->{'data'},$s);$total-=$s;
				if($total==0){last;}
			}
			if($total>0){
				# Problem sending the data -- it didn't all get sent.
				carp("respond_to_box(): problem sending data to $box.\n");
				return undef;
			}
			# We have used this request.
			$box_data->{'req'}='';
			# We are done here -- signal if we are done completely
			if(length($box_data->{'data'})==0){return 1;}
		}
		# We are done for now, but not done sending data.
		return 0;
	}else{
		carp("respond_to_box(): unknown box type ".$box_data->{'action_type'}." for $box.\n");
		return undef;
	}
}


# PARAMETER MANIPULATION ##############


# Declare a scalar to be a parameter (and defines it) which will be
# recorded with each run.  Just give an array of names (strings), and this
# will add the package names if not given, so they can be referenced later.
# Note that for this to work, the variable must be declared with "our"
# instead of "my".
# Note that duplication is allowed.
sub declare_parameter(@){
	my $pkg=caller();
	foreach my $n (@_){
		if($n=~/::/){push(@parameters,$n);}
		else{push(@parameters,"$pkg\::$n");}
	}
}


# Dump parameters (and date) to stdout, or to a file if a name is given.
# Returns true for success, and undef for error.
sub dump_parameters(;$){
	my($file)=@_;
	my($P,$n);
	# Open the write handle.
	if(!defined($file)){
		if(!open($P,'>&STDOUT')){
			carp("dump_parameters(): unable to duplicate stdout ($!).\n");
			return undef;
		}
	}else{
		if(!open($P,'>',$file)){
			carp("dump_parameters(): unable to create parameters file ($!).\n");
			return undef;
		}
	}
	# Record date.
	my $date=localtime();
	if(!(print $P ("# $date\n"))){
		carp("dump_parameters(): error writing date ($!).\n");
		close($P);return undef;
	}
	foreach my $name (@parameters){
		my $n;
		local $,;
		# Strip package name.
		($n)=($name=~/(?:::)*([^:]*)$/);
		no strict 'refs';
		if(!(print $P ("${n}=",${$name},"\n"))){
			carp("dump_parameters(): error writing parameters ($!).\n");
			close($P);return undef;
		}
	}
	if(!close($P)){
		carp("dump_parameters(): error closing parameters handle ($!).\n");
		return undef;
	}
	return 1;
}


# Given a filename in a format similar to that given by dump_parameters(),
# this loads those parameters in (and declares them).
# Returns true for success, and undef for error.
sub load_parameters($){
	my($name)=@_;
	my $pkg=caller();
	my $fd;if(!open($fd,"<",$name)){
		carp("load_parameters(): error opening file $name ($!).\n");
		return undef;
	}
	foreach my $line (<$fd>){
		chomp($line);
		my($var,$value)=($line=~/^\s*([a-zA-Z_][a-zA-Z0-9_]*)\s*=\s*(.*)$/);
		if(!defined($var) || !defined($value)){next;}
		# Declare the parameter.
		# Don't call the subroutine because that would mess up the value of
		# caller().
		if($var=~/::/){push(@parameters,$var);}
		else{push(@parameters,"$pkg\::$var");}
		# Set the parameter.
		no strict 'refs';
		${"$pkg\::$var"}=$value;
	}
	close($fd);return 1;
}


# ACTION/EVENT MANIPULATION ###########
# An action refers to something you can tell a digital or analog box to do,
# at a given time.

# For a digital box, it is a set of bits that get set hi, and another set
# of bits that get set low.

# For an analog box, you can either set a particular output to a particular
# value, or tell it to ramp from the last set point and end at the given
# time and value.  The sets are placed as close to the given time as
# possible, and the ramps are interspersed among those.

# An event is a set of actions grouped together (possibly different boxes
# and times).

# Creates a digital action -- that is, a list of bits (array ref) to set
# high and a list of bits (array ref) to set low, for a given named box.
# The arrays are copied in this routine, so modifying the passed arrays
# later won't affect these.
# The first argument is the name of the box.
# The second argument is an array reference of bits to set high (0, 1, 2,
#   etc.).
# The third argument is an array reference of bits to set low.
# An action only affects one box -- if you want to affect multiple boxes,
# use an event (a collection of actions).
# Returns the reference to the Digital_Action (undef for error).
# Will warn you if you have bad arguments.
sub digital_action($$$){
	my($box,$hi_ref,$lo_ref)=@_;
	my($bit);
	# Sanity checks and initialization.
	# All actions are defined to be at time 0 -- add them to an event to set
	# the time.
	if(ref($hi_ref) ne "ARRAY" || ref($lo_ref) ne "ARRAY"){
		carp("digital_action() needs array references.\n");return undef;
	}
	if(ref($box) ne ""){
		carp("digital_action() needs to be given a name.\n");return undef;
	}
	my $dig_ref={'time'=>0,'box'=>$box,'hi'=>[],'lo'=>[]};
	# Add the hi bits, with lots of sanity checks.
	foreach $bit (@{$hi_ref}){
		if($bit!=int($bit)){
			carp("digital_action() needs integer bit values (not $bit).\n");next;
		}
		if($bit<0 || $bit>31){
			carp("digital_action() takes bits from 0 to 31, not $bit.\n");next;
		}
		if(grep(/^$bit$/,@{$dig_ref->{'hi'}})){
			carp("digital_action(): tried to reset a hi bit hi ($bit).\n");next;
		}
		# This is not necessary, but here in case I change other stuff --
		# currently, the lo array should be empty.
		if(grep(/^$bit$/,@{$dig_ref->{'lo'}})){
			carp("digital_action(): ignoring setting bit $bit lo (already set lo).\n");
			next;
		}
		push(@{$dig_ref->{'hi'}},$bit);
	}
	# Add the lo bits, with lots of sanity checks.
	foreach $bit (@{$lo_ref}){
		if($bit!=int($bit)){
			carp("digital_action() needs integer bit values (not $bit).\n");next;
		}
		if($bit<0 || $bit>31){
			carp("digital_action() takes bits from 0 to 31, not $bit.\n");next;
		}
		if(grep(/^$bit$/,@{$dig_ref->{'lo'}})){
			carp("digital_action(): tried to reset a lo bit lo ($bit).\n");next;
		}
		if(grep(/^$bit$/,@{$dig_ref->{'hi'}})){
			carp("digital_action(): ignoring setting bit $bit lo (already set hi).\n");
			next;
		}
		push(@{$dig_ref->{'lo'}},$bit);
	}
	# Bless and return.
	bless($dig_ref,"Digital_Action");
}


# Digital_Action s have list references.  When I copy an action, I don't
# want somebody to change the original list, and affect all copies of the
# action, so this special copy routine allow makes a copy of the lists, as
# well.  Copy routines are pretty simple:  Given a hash reference, and the
# action to copy.  This will insert data in the first hash reference (not
# deleting anything already there) from the second, and bless the first.
# Returns true for success, undef for error.
sub dig_act_copy($$){
	my($dst,$src)=@_;
	if(ref($src) ne "Digital_Action"){
		carp("dig_act_copy(): Told to copy something other that a Digital_Action.\n");
		return undef;
	}
	$dst->{'time'}=$src->{'time'};
	$dst->{'box'} =$src->{'box'};
	# These makes separate copies of the arrays that 'hi' and 'lo' reference,
	# so that if the originals change, these copies won't.
	$dst->{'hi'}  =[@{$src->{'hi'}}];
	$dst->{'lo'}  =[@{$src->{'lo'}}];
	bless($dst,"Digital_Action");
}


# Creates an analog set action.  Give it a box name, a channel, and the
# voltage to set it at, and this will create an analog action corresponding
# to that.
# The first argument is the name of the box.
# The second argument is the channel number to use.
# The third argument is the voltage to set it at.
# An action here only affects one channel on one box.  Use an event to
# affect more than that.
# Returns the reference to the Analog_Action (undef for error).
# Will warn you if you have bad arguments.
# NOTE: analog outputs will not be initialized to anything special unless
# you explicitly set them.
sub analog_action_set($$$){
	my($box,$out,$val)=@_;
	# Sanity checks and initialization.
	# All actions are defined to be at time 0 -- add them to an event to set
	# the time.
	if($out!=int($out)){
		carp("analog_action_set() needs integer channel numbers (not $out).\n");
		return undef;
	}
	if($out<0 || $out>15){
		carp("analog_action_set() takes channels from 0 to 15, not $out.\n");
		return undef;
	}
	if($val<-10 || $val>10){
		carp("analog_action_set() takes values from -10 to 10, not $val.\n");
		return undef;
	}
	my $ana_ref={'time'=>0,'box'=>$box,'type'=>'set','out'=>$out,'val'=>$val};
	# Bless and return.
	bless($ana_ref,"Analog_Action");
}


# Creates an analog ramp action.  Give it a box name, a channel, and the
# voltage to end it at, and this will create an analog action corresponding
# to that.
# The first argument is the name of the box.
# The second argument is the channel number to use.
# The third argument is the voltage to end at.
# An action here only affects one channel on one box.  Use an event to
# affect more than that.
# Returns the reference to the Analog_Action (undef for error).
# Will warn you if you have bad arguments.
# NOTE: analog outputs will not be initialized to anything special unless
# you explicitly set them.
sub analog_action_ramp($$$){
	my($box,$out,$val)=@_;
	# Sanity checks and initialization.
	# All actions are defined to be at time 0 -- add them to an event to set
	# the time.
	if($out!=int($out)){
		carp("analog_action_ramp() needs integer channel numbers (not $out).\n");
		return undef;
	}
	if($out<0 || $out>15){
		carp("analog_action_ramp() takes channels from 0 to 15, not $out.\n");
		return undef;
	}
	if($val<-10 || $val>10){
		carp("analog_action_ramp() takes values from -10 to 10, not $val.\n");
		return undef;
	}
	my $ana_ref={'time'=>0,'box'=>$box,'type'=>'ramp','out'=>$out,'val'=>$val};
	# Bless and return.
	bless($ana_ref,"Analog_Action");
}


# Creates an event.  Give it a time and then a list of actions and events,
# and it returns the event reference (undef for error).  The time is
# in milliseconds.  What really happens is the given actions and events are
# expanded to just actions, and all their internal times are incremented by
# the given time.  The actions are copied, so nothing inside this event is
# connected to anything outside this event.
# Order is preserved.
# If you really want, you can have empty events (just a time).
sub event($;@){
	my $time=shift;
	my @event=();
	# First, go through the list, adding any actions to @action, and
	# expanding Events into actions.
	# Actually copy actions into @action, so that all these actions are
	# independent from other things.
	while(scalar(@_)){
		my $action=shift(@_);
		my $action_type=ref($action);
		if($action_type eq "Event"){
			unshift(@_,@{$action});
		}elsif(defined($Action{$action_type})){
			my %action;
			if(defined($Action{$action_type}{'copy'})){
				my $ret=&{$Action{$action_type}{'copy'}}(\%action,$action);
				if(!defined($ret)){
					carp("event(): could not copy data.\n");
					return undef;
				}
			}else{%action=%{$action};}
			push(@event,bless(\%action,$action_type));
		}else{
			carp("event(): skipping something that does not appear to be an action or event ($action_type).\n");
		}
	}
	# Now, go through what now should be a list of actions, and increment all
	# the times.
	foreach my $action (@event){
		$action->{'time'}+=$time;
	}
	# Bless and return.
	bless(\@event,"Event");
}


# If a name is given, sets the name of the master box to the given name,
# and returns the previous value.  Otherwise, just returns the current
# value.  Use an empty string to clear the master box (not too sure why you
# would want to, though).
sub master_box(;$){
	my($name)=@_;
	my $prev=$master_box;
	if(defined($name)){$master_box=$name;}
	return $prev;
}


# An explanation is in order of the subroutine below.  It uses a lot of
# complicated hashes, which I'll document here:
# Digital_Action (blessed as "Digital_Action"):
# The keys are:
#   'time' => the time this action takes place (in milliseconds)
#   'box'  => the name of the box to program
#   'hi'   => a list reference of the bits to set high for this action
#   'lo'   => a list reference of the bits to set lo for this action
#
# Analog_Action (blessed as "Analog_Action"):
#   'time' => the time this action takes place (in milliseconds)
#   'box'  => the name of the box to program
#   'type' => determines the type of action:
#             'set'  -- sets the output voltage to val at this time.
#             'ramp' -- generate a smooth ramp from the last set value
#                       to the given value, ending at 'time'.
#   'out'  => 0-15, which output this affects
#   'val'  => value, in volts, for the output (output is just set to this)
#   'valint' => the scaled value for val -- only used in temporary lists
#               inside certain subroutines
# Other Action types may be defined by add-ons.  They would set the
# important data in the %Action hash in the GLOBALS section.
#
# Event:
#   This is actually an array reference.  It contains a list of actions.
#
# Finally, for programming the various boxes, it is convenient to make one
# more data type, this one internal to the following subroutines and some
# of the things it calls.  It isn't officially named, but I'll call it
# box_data.  The keys are each names of boxes, and the values are hash
# references.  Each internal hash may contain the following keys:
#   'action_type' => type of box.  Valid types are the keys to %Action.
#   'data' => data to send to the box (not yet sent).
#   'sent' => data already sent to the box.
#   'sock' => a socket for connecting to the box.
#   'req'  => for storing responses from the box until they tally up to a
#             full request.
#   'init' => stuff to store to a <box>.init file.  'data' gets saved
#             to a <box>.prog file so afterwards, you can see what
#             happened.  'init', if it exists, will be written to
#             <box>.init, so you can see what things were initialized as.
#             If you have this, be sure it contains enough information to
#             completely reconstruct what the initial state of the box is.
#   'UDP'  => if defined and true, this box needs to be UDP pinged.  This
#             is ignored for the master box, which will always get a ping.
#   others specific to the type of box (like maybe an 'idata' for an easier
#   to program version of the 'init' stuff).
#   Add-ons will need to make their own versions of this, using the above
#   keys as well as any they need on their own.

# Given a list of actions and events, this figures out which boxes need to
# be programmed, computes the program data for them, and, if it is
# different than the last thing this routine programmed them with,
# reprograms them.  Attempts to write the stuff it programmed to a filename
# that's the box address appended with ".prog".  If there's an 'init' key
# in the box_data hash, then it will attempt to write the associated string
# to a file whose name is the box address appended with a ".init".  That is
# intended to carry information about how the box is to be initialized.
# By default, all boxes are assumed to start in some sort of "off" state.
# To force starting in some other position, you are allowed to give events
# happening at negative times.  The initial state at time 0, when the boxes
# start, will then be the sum total of these previous events.  If
# initialization is required to get to this state, then hopefully there's
# an initialization routine defined that will handle that.  Boxes are
# always initialized (but only reprogrammed if the program has changed).
# It then sends the UDP ping to every box with the UDP flag set (except a
# master box).
# If two conflicting things end up happening at the same time, then the one
# that occurs later in the arguments will take precedence (but a warning
# will be issued).
# Returns true if everything was successful (or close enough), and undef if
# there was a bad enough error.
sub program_boxes(@){
	my($box,$action);

	# First, expand everything into an array of actions.  I will add to the
	# actions themselves, but will not change anything below that level, so I
	# will do some minor copying.  If a copy routine is available, I'll go
	# ahead and use that, just to be safe.  I give each action an 'ord'
	# value, giving its position in the original array, so that, in the event
	# of simultaneous events, the array can be sorted so that the original
	# ordering is preserved.  Other routines will then give the last one
	# precedence.  Or, at least, they should.  Add-ons may have their own
	# ideas of what to do.
	my @action=();
	$box=0;
	while(scalar(@_)){
		$action=shift(@_);
		my $action_type=ref($action);
		if($action_type eq "Event"){
			unshift(@_,@{$action});
		}elsif(defined($Action{$action_type})){
			my %action;
			if(defined($Action{$action_type}{'copy'})){
				my $ret=&{$Action{$action_type}{'copy'}}(\%action,$action);
				if(!defined($ret)){
					carp("program_boxes(): could not copy data.\n");
					return undef;
				}
			}else{%action=%{$action};}
			$action{'ord'}=$box;++$box;
			push(@action,bless(\%action,$action_type));
		}else{
			carp("program_boxes(): skipping something that does not appear to be an action or event ($action_type).\n");
		}
	}
	# Sort the actions by time, and resolve conflicts with original ordering
	# in the argument list.
	my @action_list=sort {
		my $ret;
		$ret=($a->{'time'} <=> $b->{'time'});if($ret){return $ret;}
		return($a->{'ord'} <=> $b->{'ord'});
	} @action;

	# Make the programs, sorted by box (do it type-by-type).
	# %box_data is of the box_data type of structure talked about in the
	# comments above this routine.
	my %box_data=();
	{
		foreach my $action_type (keys(%Action)){
			if(!defined($Action{$action_type}{'prog_data'})){
				carp("program_boxes(): no program data routine for action type $action_type.\n");
				return undef;
			}
			my $box_data=&{$Action{$action_type}{'prog_data'}}(@action_list);
			if(!defined($box_data)){
				carp("program_boxes(): error handling $action_type data.\n");
				return undef;
			}
			# Concatenate into %box_data, checking for duplicate names.
			foreach $box (keys(%{$box_data})){
				if(defined($box_data{$box})){
					carp("program_boxes(): $box seems to have multiple types.\n");
					return undef;
				}
				$box_data{$box}=$box_data->{$box};
			}
		}
	}

	# Save the data to logs, and make a list of boxes that need programming.
	my @need_programming=();
	foreach $box (keys(%box_data)){
		# Only write program data and mark as something to program if its a
		# programmable box ('port' is set) and if there's program data
		# ('data'):
		if(defined($Action{$box_data{$box}{'action_type'}}{'port'}) &&
		   defined($box_data{$box}{'data'})){
			if(!defined($programs{$box}) || $programs{$box} ne $box_data{$box}{'data'}){
				# Save the programs I intend to write.
				if(!open(LOG,">","$box.prog")){
					carp("program_boxes(): error opening log for $box ($!).\n");
				}else{
					binmode(LOG);
					if(!print LOG ($box_data{$box}{'data'})){
						carp("program_boxes(): error writing log for $box ($!).\n");
					}
					if(!close(LOG)){
						carp("program_boxes(): error closing log for $box ($!).\n");
					}
				}
			}
			push(@need_programming,$box);
		}
		# Save initialization data, if there is some.
		if(defined($box_data{$box}{'init'})){
			if(!open(LOG,">","$box.init")){
				carp("program_boxes(): error opening init log for $box ($!).\n");
			}else{
				binmode(LOG);
				if(!print LOG ($box_data{$box}{'init'})){
					carp("program_boxes(): error writing init log for $box ($!).\n");
				}
				if(!close(LOG)){
					carp("program_boxes(): error closing init log for $box ($!).\n");
				}
			}
		}
	}

	# The remainder of the routine is just talking to boxes, which should be
	# skipped if we are just testing:
	if($testing){return 1;}

	# Initialize every box, if applicable.  Initializations are not done in
	# parallel -- they are assumed to go fairly fast.
	foreach $box (keys(%box_data)){
		if(defined($Action{$box_data{$box}{'action_type'}}{'init'})){
			my $ret=&{$Action{$box_data{$box}{'action_type'}}{'init'}}($box_data{$box},$box);
			if(!defined($ret)){
				carp("program_boxes(): error intializing $box.\n");
				return undef;
			}
		}
	}

	# Send data as it is requested.  Keep a list of the boxes that are
	# currently being programmed:
	my %being_programmed=();
	for(;;){
		my $data='';
		my $alive=0;
		# Open connections:
		while(scalar(keys(%being_programmed))<$max_socket &&
		      scalar(@need_programming)>0){
			$box=pop(@need_programming);
			# When the connections are opened, that clears the memory of the box,
			# so delete the appropriate value of %programs (we will set it again
			# when we have successfully written stuff -- that will happen unless I
			# encounter an error, in which case I'd rather assume I didn't know
			# what was in the box, which is why I leave the value deleted).
			if(!defined($Action{$box_data{$box}{'action_type'}}{'port'})){
				carp("program_boxes(): no port defined for $box ($box_data{$box}{'action_type'}).\n");
				foreach $box (keys(%being_programmed)){
					close($box_data{$box}{'sock'});
					delete($box_data{$box}{'sock'});
					delete($being_programmed{$box});
				}
				return undef;
			}
			$box_data{$box}{'sock'}=open_socket($box,$Action{$box_data{$box}{'action_type'}}{'port'});
			if(!defined($box_data{$box}{'sock'})){
				carp("program_boxes(): error connecting to $box.\n");
				foreach $box (keys(%being_programmed)){
					close($box_data{$box}{'sock'});
					delete($box_data{$box}{'sock'});
					delete($being_programmed{$box});
				}
				return undef;
			}
			delete($programs{$box});$being_programmed{$box}=1;
			# Some boxes may need you to send something first, so call the
			# respond routine with undef response to give it the chance to do
			# that.
			if(!defined($Action{$box_data{$box}{'action_type'}}{'respond'})){
				carp("program_boxes(): no respond routine defined for no port defined for $box ($box_data{$box}{'action_type'}).\n");
				foreach $box (keys(%being_programmed)){
					close($box_data{$box}{'sock'});
					delete($box_data{$box}{'sock'});
					delete($being_programmed{$box});
				}
				return undef;
			}
			my $ret=&{$Action{$box_data{$box}{'action_type'}}{'respond'}}($box_data{$box},undef,$box);
			if(!defined($ret)){
				carp("program_boxes(): error initializing with $box.\n");
				foreach $box (keys(%being_programmed)){
					close($box_data{$box}{'sock'});
					delete($box_data{$box}{'sock'});
					delete($being_programmed{$box});
				}
				return undef;
			}
		}
		# Form a select list and select.  A timeout is bad.  If there are no
		# sockets left, we are done.
		if(scalar(keys(%being_programmed))==0){last;}
		my $vec='';foreach $box (keys(%being_programmed)){
			vec($vec,fileno($box_data{$box}{'sock'}),1)=1;
		}
		select($vec,undef,undef,$sock_timeout);
		# Read everything that wants reading and handle requests from boxes
		# (closing sockets when done on either end -- warn if the other end
		# closed first).
		$alive=0;foreach $box (keys(%being_programmed)){
			# Skip boxes that are not ready to read.
			if(!vec($vec,fileno($box_data{$box}{'sock'}),1)){next;}
			++$alive;
			sysread($box_data{$box}{'sock'},$data,65536);
			if(length($data)==0){
				#	This seems to be a sign that the connection was closed.
				carp("program_boxes(): the connection to $box was unexpectedly closed.\n");
				close($box_data{$box}{'sock'});
				delete($box_data{$box}{'sock'});
				delete($being_programmed{$box});
			}
			my $ret=&{$Action{$box_data{$box}{'action_type'}}{'respond'}}($box_data{$box},$data,$box);
			if(!defined($ret)){
				carp("program_boxes(): error responding to $box.\n");
				foreach $box (keys(%being_programmed)){
					close($box_data{$box}{'sock'});
					delete($box_data{$box}{'sock'});
					delete($being_programmed{$box});
				}
				return undef;
			}elsif($ret<0){
				carp("program_boxes(): $box filled up.\n");
				close($box_data{$box}{'sock'});
				delete($box_data{$box}{'sock'});
				delete($being_programmed{$box});
			}elsif($ret>0){
				# Successfully done.  Save what we sent to the box.
				close($box_data{$box}{'sock'});
				delete($box_data{$box}{'sock'});
				delete($being_programmed{$box});
				$programs{$box}=$box_data{$box}{'sent'};
			}
			# Otherwise, we are still going with this box.
		# End of reading any available data from boxes.
		}
		if(!$alive){
			my $strg="program_boxes(): timed out waiting for:\n";
			foreach $box (keys(%being_programmed)){
				close($box_data{$box}{'sock'});
				delete($box_data{$box}{'sock'});
				delete($being_programmed{$box});
				$strg.="$box\n";
			}
			carp($strg);return undef;
		}
	# End of sending data to boxes.
	}

	# Give the boxes a little time to process their data.
	Time::HiRes::sleep(0.01);
	# Send the UDP pings.  Do not ping the master box, if there is one.
	{
		my %udp_box_data=%box_data;
		if(defined($master_box) && defined($udp_box_data{$master_box})){
			delete($udp_box_data{$master_box});
		}
		# Don't ping boxes without a true 'UDP' value:
		foreach $box (keys(%udp_box_data)){
			if(!$udp_box_data{$box}{'UDP'}){
				delete($udp_box_data{$box});
			}
		}
		my $udp_rsp=udp_ping(\%udp_box_data);
		my $err=undef;
		if(!defined($udp_rsp)){
			carp("program_boxes(): unable to send UDP ping.\n");
			return undef;
		}
		foreach $box (keys(%udp_box_data)){
			if(!defined($udp_rsp->{$box})){
				$err.="program_boxes(): $box did not respond to UDP ping.\n";
			}elsif(!$udp_rsp->{$box}){
				$err.="program_boxes(): $box gave a bad response to UDP ping\n";
			}
		}
		if(defined($err)){
			carp($err);return undef;
		}
	}

	# We are finally done.  Store which boxes were programmed, and return.
	@programmed_boxes=keys(%box_data);
	# Wait a tiny amount, so the boxes can get ready for a UDP ping.
	Time::HiRes::sleep(0.001);
	return 1;
}


# Sends the "go" signal to master box, which should, hopefully, trigger the
# other boxes.
# Returns true for success, and undef for error.
sub trigger(){
	if(!defined($master_box)){
		carp("trigger(): master_box name is not defined.\n");
		return undef;
	}
	if($master_box eq ''){
		carp("trigger(): master_box name is blank.\n");
		return undef;
	}

	# If we are just testing, return straight away.
	if($testing){return 1;}

	# Two cases -- if we just programmed it, then send it a ping and assume
	# it will trigger whatever is necessary.  Otherwise, connect through the
	# control port and set a channel high.
	my $programmed=0;
	foreach my $box (@programmed_boxes){
		if($box eq $master_box){$programmed=1;last;}
	}
	if($programmed){
		my $udp_rsp=udp_ping({$master_box=>undef});
		my $err=undef;
		if(!defined($udp_rsp)){
			carp("trigger(): unable to send UDP ping to $master_box.\n");
			return undef;
		}
		if(!defined($udp_rsp->{$master_box})){
			carp("trigger(): $master_box did not respond to UDP ping.\n");
			return undef;
		}elsif(!$udp_rsp->{$master_box}){
			carp("trigger(): $master_box gave a bad response to UDP ping\n");
			return undef;
		}
		return 1;
	}else{
		my $sock=open_socket($master_box,$control_port);
		if(!defined($sock)){
			carp("trigger(): unable to open socket to $master_box.\n");
			return undef;
		}
		if(!expect($sock,"","200") ||  # Get the welcome message.
		   !expect($sock,"o","400") ||
		   !expect($sock,"V","412") ||
		   !expect($sock,"1","Voltage") ||
		   !expect($sock,"5","set to") ||
		   !expect($sock,"V","412") ||
		   !expect($sock,"1","Voltage") ||
		   !expect($sock,"0","set to") ||
		   !expect($sock,"q","200")){
			carp("trigger(): communication error with $master_box.\n");
			close($sock);return undef;
		}
		write_socket($sock,"q\r\n");close($sock);
		return 1;
	}
}


# Sends one ping to each of the boxes that was most recently programmed,
# and waits for a response.  The optional argument overrides the default
# period to wait (in seconds).
# Returns true for success, and undef for some error (like a timeout).
# Returns true automatically (no wait) if you are just testing.
sub wait_for_boxes(;$){
	my($wait)=@_;if(!defined($wait)){$wait=$sock_timeout;}
	my(%box,$box);
	my $msg='';

	if($testing){return 1;}

	# Ping all the boxes -- need to use actual ping program because that's
	# SUID.
	my $timeout=POSIX::ceil($wait);
	foreach $box (@programmed_boxes){
		my $command={
			# -c 1 means only one packet is sent.
			# -t <timeout> means to quit after waiting timeout seconds.
			# -i <delay> means to wait delay between packets.  For some reason,
			# when -t is large, and -c is 1, it won't wait long enough without a
			# -i option, too.
			'command'=>[$ping_prog,'-c','1','-i',$timeout,'-t',$timeout,$box],
			'indirect'=>$ping_prog,
			'retval'=>-1,
			'pid'=>-1,
		};
		if(!($ping_fork->queue($command))){
			$msg.="  $ping_fork->{'error'}\n";
			last;
		}
		$box{$box}=$command;
	}
	# Now, they've all been pinged.  Wait a bit, and check for return values.
	# Need to wait a bit longer than timeout so that ping can return.
	# Do that by sleeping for $timeout first, and then keep sleeping that
	# interval until processes stop finishing (sleep should stop when the
	# first child dies, so this doesn't wait overly long).  Give it a little
	# extra (double) since ping isn't too good about returning on time.
	sleep($timeout);
	{
		my @status=$ping_fork->status();
		my $left=$status[0]+$status[1];
		my $time=Time::HiRes::time();
		while($left>0){
			sleep($timeout);
			@status=$ping_fork->status();
			if($left!=$status[0]+$status[1]){
				$left=$status[0]+$status[1];
				$time=Time::HiRes::time();
			}
			if(Time::HiRes::time()-$time>2*$timeout){last;}
		}
	}
	# Check what's left.
	foreach $box (@programmed_boxes){
		if($box{$box}{'retval'}<0){
			# It didn't return -- shouldn't happen -- it shouldn't take this
			# long:
			$msg.="  process for $box did not quit.\n";
		}elsif($box{$box}{'retval'}!=0){
			# Returned with an error.
			my($ret,$sig)=($box{$box}{'retval'}>>8,$box{$box}{'retval'}&0xFF);
			$msg.="  $box did not respond (ping returned $ret, signal $sig).\n";
		}
		# Do nothing for a successful return (0).
	}
	# Kill anything that is left -- they are just ping programs, so do a hard
	# kill.
	$ping_fork->killall('KILL');
	if($msg eq ''){return 1;}
	carp("wait_for_boxes(): errors occurred:\n$msg");return undef;
}


# CAMERA SERVER ROUTINES ##############


# Connects to my camera server.
# Returns undef for error, true for success.
sub cam_open(){
	if(defined($cam_sock)){
		carp("cam_open(): we seem to be connected already.\n");
		return undef;
	}
	if(!socket($cam_sock,AF_UNIX,SOCK_STREAM,0)){
		carp("cam_open(): unable to create socket ($!).\n");
		$cam_sock=undef;return undef;
	}
	if(!connect($cam_sock,sockaddr_un($cam_socket))){
		carp("cam_open(): unable to connect to socket ($!).\n");
		carp("            If camserver is not running, run it.\n");
		carp("            If it won't start because the socket exists,\n");
		carp("            delete the socket and start camserver.\n");
		carp("            Disconnecting the USB or cycling the power on the\n");
		carp("            camera tend to disconnect the USB, and the Mac USB\n");
		carp("            toolkit I use won't reconnect, so kill -15 the\n");
		carp("            camserver and restart it.\n");
		close($cam_sock);$cam_sock=undef;return undef;
	}
	{my $oldfh=select($cam_sock);binmode($cam_sock);$|=1;select($oldfh);}
	# Get the PID.
	if(!cam_command("PID")){
		carp("cam_open(): error sending PID command.\n");
		return undef;
	}
	# Wait for a response, but don't wait too long.
	my $end_time=Time::HiRes::time()+$sock_timeout;
	my($bits,$data,$total);$total='';$cam_PID=undef;
	while((my $time=Time::HiRes::time())<$end_time && !defined($cam_PID)){
		$bits='';vec($bits,fileno($cam_sock),1)=1;
		select($bits,undef,undef,$end_time-$time);
		if(vec($bits,fileno($cam_sock),1)){
			sysread($cam_sock,$data,65536);
			$total.=$data;my(@rsp)=split(/\n+/,$total);
			# Get the PID:
			foreach(@rsp){(/^(\d+)$/) && ($cam_PID=$1);}
			# Only keep the last line:
			$total=~s/^.*\n//os;
		}
	}
	return 1;
}


# Sends a string to my camera server (appends the '\n' for you).
# In the event of an error, returns undef and the camera socket will be
# closed.  Returns true for success.
sub cam_command($){
	my($com)=@_;
	if(!defined($cam_sock)){
		carp("cam_command(): camera socket does not seem to be opened.\n");
		return undef;
	}
	dump_handle($cam_sock);
	if(!defined(syswrite($cam_sock,"$com\n"))){
		carp("cam_command(): error sending to socket.  Closing socket.\n");
		close($cam_sock);$cam_sock=undef;return undef;
	}
	return 1;
}


# Sends the takepic command to the camera -- requires only the filename.
# Prepends the current directory if the name doesn't start with a '/'.
# Returns true for success, negative for error (waits a bit for a response
# from the camera that it is waiting).
# If we are in testing mode, this just returns (after checking that the
# connection is open).
sub cam_takepic($){
	my($path)=@_;
	if(!defined($cam_sock)){
		carp("cam_takepic(): camera socket does not seem to be opened.\n");
		return undef;
	}
	dump_handle($cam_sock);
	if($testing){return 1;}
	if(!($path=~/^\//)){
		my $pwd=`$pwd_prog`;chomp($pwd);
		$path="${pwd}/${path}";
	}
	if(!cam_command("takepic $path")){
		carp("cam_takepic(): error sending picture command.\n");
		return undef;
	}
	# Wait for a response, but don't wait too long.
	my $end_time=Time::HiRes::time()+$sock_timeout;
	my($bits,$data,$total);$total='';
	while((my $time=Time::HiRes::time())<$end_time){
		$bits='';vec($bits,fileno($cam_sock),1)=1;
		select($bits,undef,undef,$end_time-$time);
		if(vec($bits,fileno($cam_sock),1)){
			sysread($cam_sock,$data,65536);
			$total.=$data;my(@rsp)=split(/\n+/,$total);
			# Here's the response I wanted:
			if(grep(/^Taking picture/,@rsp)){return 1;}
			# Only keep the last line:
			$total=~s/^.*\n//os;
		}
	}
	carp("cam_takepic(): timed out waiting for picture response.\n");
	return undef;
}


# Waits for the camera server to respond again (may time out).  If you
# always use this sometime before calling another cam_* command after using
# cam_takepic(), you should never get out of sync with the camera server.
# The optional argument overrides the default timeout (in seconds).
# In the event of a timeout, this tries to cancel any pending exposures, so
# you can just quit the program that called this, or not, as you'd like,
# but further attempts at communication with the server should succeed.
# If you are in testing mode, this just returns after checking that the
# connection is open.
sub cam_wait(;$){
	my($wait)=@_;if(!defined($wait)){$wait=$sock_timeout;}
	if(!defined($cam_sock)){
		carp("cam_wait(): camera socket does not seem to be opened.\n");
		return undef;
	}
	dump_handle($cam_sock);
	if($testing){return 1;}
	# Use the hello command for this.
	if(!cam_command("hello")){
		carp("cam_wait(): error sending command.\n");
		return undef;
	}
	# Wait for a response, but don't wait too long.
	my $end_time=Time::HiRes::time()+$wait;
	my($bits,$data,$total);$total='';
	while((my $time=Time::HiRes::time())<$end_time){
		$bits='';vec($bits,fileno($cam_sock),1)=1;
		select($bits,undef,undef,$end_time-$time);
		if(vec($bits,fileno($cam_sock),1)){
			sysread($cam_sock,$data,65536);
			$total.=$data;my(@rsp)=split(/\n+/,$total);
			# Here's the response I wanted:
			if(grep(/^Hello/,@rsp)){return 1;}
			# Only keep the last line:
			$total=~s/^.*\n//os;
		}
	}
	carp("cam_wait(): timed out waiting for picture response.\nWill try to cancel pending exposures.\n");
	if(defined($cam_PID)){kill('INT',$cam_PID);}
	return undef;
}


# Disconnects from my camera server (and cancels pending exposures).
# Returns true.
sub cam_close(){
	if(defined($cam_PID)){kill('INT',$cam_PID);}
	if(defined($cam_sock)){close($cam_sock);}
	$cam_sock=undef;$cam_PID=undef;return 1;
}


# DISPLAY ROUTINES ####################


# Displays an image using the image server.
# The only required argument is the filename to display (if you make it
# undef, this will try to start a server and set it up with whatever
# size and position you give).
# The optional argument is a hash that sets various options (all have
# defaults or aren't used if not given).  Keys are:
#   false => true or false, depending on whether you want false colors.
#   pos   => 2-element array reference, with x-y position of window.
#   size  => 2-element array reference, with widht-height of window
#            (negative values should match image size).
#   name  => Name of window (so you can have multiple windows).
# This works by running the image server in the background with the
# Fork_System library.
# Returns true for success (or so it thinks) and undef for an error
# (probably not worth quitting over, though -- just warn).
sub display($;%){
	my($path,%options)=@_;
	# Make an absolute path, if warranted.
	if(defined($path) && !($path=~/^\//)){
		my $pwd=`$pwd_prog`;chomp($pwd);
		$path="${pwd}/${path}";
	}
	# Set defaults.
	if(!defined($options{'false'})){$options{'false'}=0;}
	my $name;
	if(!defined($options{'name'})){$name=$default_img_name;}
	else{$name=$options{'name'};}
	# If the program does not seem to be open, or seems to have died, run a
	# new program.
	if(!defined($img_comm{$name}) || $img_comm{$name}{'retval'}>=0){
		if(defined($img_comm{$name}{'stdin'})){close($img_comm{$name}{'stdin'});}
		if(defined($img_comm{$name}{'stdout'})){close($img_comm{$name}{'stdout'});}
		my $command={
			'command'=>[$img_prog,$name],
			'indirect'=>$img_prog,
			'retval'=>-1,
			'pid'=>-1,
			'stdin'=>'',
			'stdout'=>'',
			# Ignore stderr.
		};
		if(!($img_fork->queue($command))){
			carp("display(): error: $img_fork->{'error'}.\n");
			return undef;
		}
		$img_comm{$name}=$command;
	}
	# Now, send the commands.
	if(!dump_handle($img_comm{$name}{'stdout'})){
		carp("display(): unable to read from pipe.\n");
		close($img_comm{$name}{'stdin'});
		close($img_comm{$name}{'stdout'});
		delete($img_comm{$name});
		return undef;
	}
	my $command='';
	if(defined($path)){
		if(defined($options{'false'}) && $options{'false'}){
			$command="false $path\n";
		}else{$command="display $path\n";}
	}
	if(ref($options{'pos'}) eq 'ARRAY' && scalar(@{$options{'pos'}})==2){
		$command.="move $options{'pos'}[0] $options{'pos'}[1]\n";
	}
	if(ref($options{'size'}) eq 'ARRAY' && scalar(@{$options{'size'}})==2){
		$command.="resize $options{'size'}[0] $options{'size'}[1]\n";
	}
	# Don't bother selecting on the write -- if the pipe is closed, the write
	# will fail (but do block the SIGPIPE). 
	# Also, don't bother waiting for the response -- let it display on its
	# own while we go on.
	local $SIG{'PIPE'}='IGNORE';
	if(!syswrite($img_comm{$name}{'stdin'},$command)){
		carp("display(): error sending command to image server ($!).\n");
		close($img_comm{$name}{'stdin'});
		close($img_comm{$name}{'stdout'});
		delete($img_comm{$name});
		return undef;
	}
	dump_handle($img_comm{$name}{'stdout'});
	return 1;
}


# Closes the window and image server with the given name (same name used
# for display above).
# Returns true for success, 0 if there's nothing by that name, and undef
# for error (probably not a bad error, though).
sub close_display($){
	my($name)=@_;
	if(defined($img_comm{$name})){
		# It should be enough to just close stdin and send a nice kill to shut
		# down the image server.
		if(defined($img_comm{$name}{'stdin'})){close($img_comm{$name}{'stdin'});}
		if(defined($img_comm{$name}{'stdout'})){close($img_comm{$name}{'stdout'});}
		if(defined($img_comm{$name}{'pid'}) && $img_comm{$name}{'pid'}>0){
			kill('TERM',$img_comm{$name}{'pid'});
		}
		delete($img_comm{$name});
		return 1;
	}
	return 0;
}


# COMMAND ROUTINES ####################
# These are meant to take the place of the system call.


# This runs a command in the background (given the command and the
# arguments).  The stdout and stderr are left alone.
# Returns true for success, and undef for error.
sub queue_com(@){
	my $command={
		'command'=>\@_,
		'pid'=>-1,
		'retval'=>-1,
		'stdout'=>undef,
		'stderr'=>undef,
	};
	if(!defined($user_fork->queue($command))){
		carp("queue_com(): error queuing process: $user_fork->{'error'}.\n");
		return undef;
	}
	return 1;
}
# Same as above, but stdout and stderr are piped to /dev/null.
# Returns true for success, and undef for error.
sub queue_com_silent(@){
	my $command={
		'command'=>\@_,
		'pid'=>-1,
		'retval'=>-1,
	};
	if(!defined($user_fork->queue($command))){
		carp("queue_com_silent(): error queuing process: $user_fork->{'error'}.\n");
		return undef;
	}
	return 1;
}


# Waits for alll the commands to finish.
# Returns true for success, and undef for error.
sub wait_com(){
	while($user_fork->running()){sleep($sock_timeout);}
	return 1;
}


# Runs the given command (just like queue_com()), but waits for completion,
# and returns what the command returned.  This is intended to be like a
# system() call.
sub run_com(@){
	my $command={
		'command'=>\@_,
		'pid'=>-1,
		'retval'=>-1,
		'stdout'=>undef,
		'stderr'=>undef,
	};
	if(!defined($user_fork->queue($command))){
		carp("run_com(): error queuing process: $user_fork->{'error'}.\n");
		return undef;
	}
	while($command->{'retval'}<0){sleep($sock_timeout);}
	return $command->{'retval'};
}
# Same as above, but sends stdout and stderr to /dev/null.
sub run_com_silent(@){
	my $command={
		'command'=>\@_,
		'pid'=>-1,
		'retval'=>-1,
	};
	if(!defined($user_fork->queue($command))){
		carp("run_com_silent(): error queuing process: $user_fork->{'error'}.\n");
		return undef;
	}
	while($command->{'retval'}<0){sleep($sock_timeout);}
	return $command->{'retval'};
}
# Same as above, but returns the stdout output of the command.
sub run_com_backtick(@){
	my $command={
		'command'=>\@_,
		'pid'=>-1,
		'retval'=>-1,
		'stdout'=>1,
		'stderr'=>undef,
	};
	if(!defined($user_fork->queue($command))){
		carp("run_com_backtick(): error queuing process: $user_fork->{'error'}.\n");
		return undef;
	}
	my $stdout='';
	for(;;){
		my $bits='';vec($bits,fileno($command->{'stdout'}),1)=1;
		select($bits,undef,undef,$sock_timeout);
		if(vec($bits,fileno($command->{'stdout'}),1)){
			my($d,$r);
			$r=read($command->{'stdout'},$d,65536);
			# The man page says undefined means an error, but sometimes read
			# seems to return undef without setting $!.  In that case, the data
			# still seemed good, so I check $! too.  An undef return value and a
			# empty error are treated as the end of output.
			if(!defined($r) && $!){
				carp("run_com_backtick(): error reading from command's stdout ($!).\n");
				return undef;
			}
			# End of output:
			if(!$r){last;}
			$stdout.=$d;
		}
	}
	return $stdout;
}


# RUN ROUTINES ########################


# Initializes a run by making a directory just for the run, chdir into it,
# and recording all the parameters into the PARAMETERS file in it.
# The optional argument lets you pick a directory name to make.
# If you have a second argument, that will be the filename to write the
# parameters to.
# Returns true for success, undef for error.
sub init_run(;$$){
	my($name,$param)=@_;
	if(!defined($param)){$param=$PARAMETERS;}
	if(defined($basedir)){
		carp("init_run(): already in a run.\n");
		return undef;
	}
	my $run=$run_num+1;
	# Make the dir, and chdir.
	if(!defined($name)){$name="run.".sprintf('%04u',$run);}
	if(!mkdir($name)){
		carp("init_run(): unable to create $name ($!).\n");
		return undef;
	}
	my $dir=`$pwd_prog`;chomp($dir);if(!defined($dir)){
		carp("init_run(): unable to get current directory.\n");
		return undef;
	}
	if(!chdir($name)){
		carp("init_run(): unable to chdir to $name ($!).\n");
		return undef;
	}
	# Record date and parameters.
	if(!dump_parameters($param)){
		carp("init_run(): unable to write parameters file ($!).\n");
		chdir($dir);return undef;
	}
	# Success.
	$basedir=$dir;$run_num=$run;return 1;
}


# Ends a run by moving out of the run directory.
# Returns true for success, undef for error.
sub end_run(){
	if(!defined($basedir)){
		carp("end_run(): no run currently initialized.\n");
		return undef;
	}
	if(!chdir($basedir)){
		carp("end_run(): unable to chdir ($!).\n");
		return undef;
	}
	$basedir=undef;
	return 1;
}


# EXPORT STUFF ########################


# If you export only subroutines, it is faster to leave off the '&' on each
# one (from the Exporter documentation).  If no type is given, the '&' is
# assumed.
our @EXPORT=qw(declare_parameter dump_parameters load_parameters
               digital_action analog_action_set analog_action_ramp event
               master_box program_boxes wait_for_boxes trigger
               cam_open cam_command cam_takepic cam_wait cam_close
               display close_display
               queue_com queue_com_silent wait_com run_com run_com_silent
               run_com_backtick
               init_run end_run);
our @EXPORT_OK=qw(reap);


}
1;
