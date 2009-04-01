#!/usr/bin/perl -w

# Author: JJT
# This provides my interface to our Agilent scope.
# It contains no built-in commands -- just routines the establish and break
# a connection, and send data back and forth (performing some parsing in
# the process).

# This is the one really intended for port 5025 and binary data transfer.
# There are no prompts, so you don't always know if the scope understood
# your command.

{ package Agilent_Scope;
use strict;
use warnings;
use Low_Socket;


# METHODS #############################


# Constructor:
#   open -- the constructor.  Give it the host and port (port optional --
#           there is a default).
#           Call it directly (Agilent_Scope::open(), not
#           Agilent_Scope->open()).

# Methods:
#   verbose -- Pass it a boolean (optional).  Returns the current verbose
#              settings.  If a value is passed, this returns the setting
#              before changing it.  No method writes to stdout or stderr if
#              verbose is false, but some write out status messages
#              (sometimes mirrored in the status string) if verbose is
#              true.
#   status -- Returns the status string.
#   send_com -- Sends the given command (not terminated).
#               Returns the scope response, if there is one.  An optional
#               second argument determines how long to wait for a response
#               (default is 0, so don't expect a response).
#   send_bin -- Sends a command that is expected to give a binary response.
#               Returns the parsed response (an array).
#   get_data -- Pass a hash, followed by a boolean and a number.
#               This returns 0 for success and undef for error.
#               See the comments for the method itself for a description
#               of the arguments.
#   close -- closes the connection.
#   DESTROY -- the destructor.  Just calls close above.

# Subroutines:
#   save_data -- Give it a filename, the channel data (ref, from get_data),
#                and the output format (0 for raw data, nonzero for
#                processed vaules), and it will save the data to a file.
#                Not a method, and not exported, so call as:
#                Agilent_Scope::save_data(...);
#                Returns positive for success, 0 or undef for error (undef
#                for real error, 0 for close error), and prints a message
#                for an error.  Also returns 0 if you give bad/empty data.

# The data used by this module is stored in a hash referenced by the class
# instance.  The keys are:
#   socket -- a reference to the Low_Socket instance used by this instance.
#   status -- a string containing the most recent status messages (from the
#             last method called).
#   verbose -- a boolean.  If true, then some of the methods may write to
#              stdout.
# Every instance of this class is assumed to be an open connection.  The
# connection is closed when the destructor is called (or the close method).
# Once closed (from either end), the methods should still return in a
# timely fashion, but may return errors.  An instance cannot be reopened --
# you have to create a new one.
# NOTE:  The methods check for some internal consistency of the internal
# data, but it is possible to cause errors by changing the internal data
# from outside those methods, so do that at your own risk.


# Send a generic command to the scope.  The command should not be
# terminated -- the '\n' is added in here.
# Returns the entire response up to a prompt (but without the prompt).
# Returns undef for some kind of error (prints an error message).
# The second argument, if included, overrides the default waiting time for
# a response (default is 0, so don't expect much response unless there's
# one from a previous command).
sub send_com($;$){
	my($this,$cmd,$timeout)=@_;

	# This needs to be a connected instance.
	if(!defined(ref($this))){return undef;}

	# Check if the connection seems to be open.
	if(!defined($this->{"socket"})){
		if($this->{"verbose"}){print "Attempted to write to unopen scope connection.\n";}
		$this->{"status"}="Fail: scope_com(): connection not open.";
		return undef;
	}

	if(!defined($timeout)){$timeout=0;}
	my $resp='';my $r;
	if(!defined($this->{"socket"}->write("$cmd\n"))){
		if($this->{"verbose"}){print "Error writing the command.\n";}
		$this->{"status"}="Failed send_com(): error writing command.";
		return undef;
	}
	$r=$this->{"socket"}->read($timeout);
	if(!defined($r)){
		if($this->{"verbose"}){print "Error ocurred while reading a response.\n";}
		$this->{"status"}="Failed send_com(): error reading response (connection closed?).";
		return undef;
	}
	return $r;
}


# This version of send_com() is specifically for the :waveform:data?
# query, although it will work for other queries that return something
# starting with the #8 stuff.  As opposed to waiting for a prompt, it
# expects a #nDDDDDDDD (where n gives the number of D's, and the D's,
# concatenated, tell how many bytes you may expect after that) and reads
# that many bytes into a string that it returns.
# Returns undef for an error.
sub send_bin($;$){
	my($this,$cmd,$timeout)=@_;

	# This needs to be a connected instance.
	if(!defined(ref($this))){return undef;}

	# Check if the connection seems to be open.
	if(!defined($this->{"socket"})){
		if($this->{"verbose"}){print "Attempted to write to unopen scope connection.\n";}
		$this->{"status"}="Fail: scope_bin(): connection not open.";
		return undef;
	}

	if(!defined($timeout)){$timeout=10;}

	my $data='';my($d,$n,$num);
	if(!defined($this->{"socket"}->write("$cmd\n"))){
		if($this->{"verbose"}){print "Error writing the command.\n";}
		$this->{"status"}="Fail: send_bin(): error writing command.";
		return undef;
	}
	$n=-1;$num=-1;for(;;){
		$d=$this->{"socket"}->read($timeout);
		if(!defined($d)){
			if($this->{"verbose"}){print "Error ocurred while reading a response.\n";}
			$this->{"status"}="Fail: send_bin(): error reading response (connection closed?).";
			return undef;
		}
		if($d eq ''){
			if($this->{"verbose"}){print "Timeout -- something's wrong.\n";}
			$this->{"status"}="Fail: send_bin(): timeout.";
			return undef;
		}
		$data.=$d;
		if($n<0 && length($data)>=length("#.")){
			# We haven't read the number of digits yet.  Look for it.
			# It may not be at the very beginning.
			# On port 5024, the scope would echo commands, but that does not seem
			# to be the case on 5025.
			if($data=~/#(\d)/){
				# Save the number.
				$n=$1;
				# Now that we do not need it, trim it off.
				$data=~s/^.*?#\d//s;
			}elsif(length($data)>1024){
				# There seems to be a problem -- too much data, nothing we want.
				if($this->{"verbose"}){print "Cannot find #n string in response.\n";}
				$this->{"status"}="Fail: send_bin(): could not find #n in response.";
				return undef;
			}
		}
		if($n>=0 && $num<0 && length($data)>=$n){
			# Grab the digits, if we can.
			$num=substr($data,0,$n);if(!($num=~/^\d*$/)){
				if($this->{"verbose"}){print "Unexpected response number.\n";}
				$this->{"status"}="Fail: send_bin(): unable to get return data size.";
				return undef;
			}
			# Trim it off.
			substr($data,0,$n)="";
		}
		# Check if we have gotten all the data we need.  If so, we are done.
		if($num>0 && length($data)>=$num){
			last;
		}
	}
	$this->{"status"}="Success: send_bin().";
	# Trim and return.
	return substr($data,0,$num);
}


# This subroutine sets up the scope to get the data in the form I usually
# want it, and then gets the data.
# The first argument is a hash ref.  The keys are the names of the channels
# you want to get data from (if the values are undef, then they are
# ignored).  The values, if defined, will be replaced with references to
# hashes.  The keys of those hashes are:
#   tinc, torig, tref, yinc, yorig, yref, data.
#   t is time, y is data.  Actual values for both are
#   (number-ref)*inc+orig.
#   For t, number is which data element.  For y, number is in the array
#   referenced by the key "data".
#   tref is supposedly always 0, and yref is supposedly steps/2, where
#   steps is 64k for format=WORD and 256 for format=BYTE.
# The second argument is a boolean for whether to digitize the requested
# channels.  This will take a new dataset.  If you choose not to digitize a
# new dataset, and the requested channels are not on, you may get an error.
# The third argument is the number of points to request.  Use "max" for
# all.
# This returns 1 for success and undef for error.
sub get_data($$$){
	my($this,$chanref,$digitize,$numpoints)=@_;

	# This needs to be a connected instance.
	if(!defined(ref($this))){return undef;}

	# This is required for :waveform commands (which I need to read out
	# data):
	defined($this->send_com(":timebase:mode main")) || return undef;
	# No averaging:
	defined($this->send_com(":acquire:type normal")) || return undef;
	# Digitize, if requested.
	if($digitize){
		my $strg=":digitize";my $comma=" ";
		foreach(keys(%{$chanref})){
			if(defined($chanref->{$_})){$strg.="$comma$_";$comma=",";}
		}
		defined($this->send_com($strg)) || return undef;
	}
	# After a :digitize, this command will not return anything until the
	# digitize is finished.  Wait a little longer than the longest scan our
	# scope has.  This should return almost immediately if there's nothing
	# going on.  Unfortunately, it also returns when you are currently
	# running a :single.
	defined($this->send_com('*OPC?',510)) || return undef;

	# Commands that are somewhat related to digitize:
	# digitize tends to stop the scope.  You can manually stop it with:
	# :stop.  You can run (:run) or do a single pass (:single), too.
	# Here, we need to make sure the scope is not waiting for another run, or
	# it may not respond to the :waveform:preamble? command.
	defined($this->send_com(':stop')) || return undef;

	foreach my $chan (keys(%{$chanref})){ 
		defined($chanref->{$chan}) || next;
		# dataref points to a new, empty hash for each channel.
		my $dataref={};
		# Set the channel from which I want information.
		defined($this->send_com(":waveform:source ".$chan)) || return undef;
		# Set number of points.
		defined($this->send_com(":waveform:points $numpoints")) || return undef;
		# I wanted 2 bytes per data point, with least-significant first, unsigned.
		#defined($this->send_com(":waveform:format word")) || return undef;
		#defined($this->send_com(":waveform:byteorder lsbfirst")) || return undef;
		# but then I found out that when you aren't averaging, there are only
		# 8-bits worth of data, so just do this:
		defined($this->send_com(":waveform:format byte")) || return undef;
		defined($this->send_com(":waveform:unsigned on")) || return undef;
		# Get and interpret the preamble.
		my $pre=$this->send_com(":waveform:preamble?",10);
		defined($pre) || return undef;
		# Remove newlines, etc. and split and assign.
		# Skip the format, type, points, and count parts -- I assume they got
		# set.
		$pre=~s/[\r\n]//g;
		if(!((undef,undef,undef,undef,$dataref->{"tinc"},$dataref->{"torig"},
		      $dataref->{"tref"},$dataref->{"yinc"},$dataref->{"yorig"},
					$dataref->{"yref"})=split(/,/,$pre))){
			if($this->{"verbose"}){print "Error parsing preamble.\n";}
			$this->{"status"}="Fail: get_data(): error parsing preamble for $chan.";
			return undef;
		}
		# Get the data.
		my $data=$this->send_bin(":waveform:data?");
		defined($data) || return undef;
		# Make an array, and place it in the hash.
		$dataref->{"data"}=[unpack("C*",$data)];
		# Place this pointer in the appropriate spot.
		$chanref->{$chan}=$dataref;
	}
	$this->{"status"}="Success: get_data().";
	return 1;
}


# The constructor.  Given a host and a port, this opens a connection, and,
# if successful, creates and returns an instance of this class.  Returns
# undef for error.  An optional third argument allows for overriding the
# default verbose setting (this routine may or may not produce output).
# Note that, since if this fails, you get nothing back, if verbose is off,
# you get nothing as to why it failed.  If really necessary, you can have
# verbose on, and trap stdout to get the error message.
sub open($$;$){
	my($host,$port,$verbose)=@_;
	my(%self,$this,$socket,$status);
	defined($verbose) || ($verbose=1);
	local $|=1;

	# Make sure this was not called like Agilent_Scope->open():
	ref($host) && return undef;

	if($verbose){print "Attempting to open a connection ... ";}
	if(!defined($socket=Low_Socket::open($host,$port))){
		if($verbose){print "Failed.\n";}return undef;
	}

	# Go ahead and form the instance, since I use it for the send_com() call.
	%self=("socket",$socket,"status",$status,"verbose",$verbose);
	$this=bless(\%self);

	if($verbose){print "Success.\n";}
	$self{"status"}="Success: open(): opened socket.";

	return $this;
}


# Sets the verbose boolean, if given, and returns the current value (before
# changes).
sub verbose(;$){
	my($this,$verbose)=@_;

	# This needs to be a connected instance.
	if(!defined(ref($this))){return undef;}

	my $ret=$this->{"verbose"};

	if(defined($verbose)){$this->{"verbose"}=$verbose;}
	return $ret;
}


# Returns the status string.
sub status(){
	my($this)=@_;

	# This needs to be a connected instance.
	if(!defined(ref($this))){return undef;}
	return $this->{"status"};
}


# Closes the connection, and undefines it internally so the routines can
# catch that.  Doesn't return anything.
sub close(){
	my($this)=@_;
	if(defined(ref($this)) && defined($this->{"socket"})){
		$this->{"socket"}->close();$this->{"socket"}=undef;
	}
}


# The destructor.  Just calls close.
sub DESTROY(){
	# Use the '&' to avoid confusion with the CORE::close() call.
	&close(@_);
}


# A function for saving data.  Currently, there are two output formats,
# both ASCII, tab-delimited formats.  The headers are basically the same,
# but one saves the original numbers ($out_format=0) and one saves the
# processed numbers ($out_format!=0).
# Give this a filename, the channel data (a reference to a hash -- same as
# for get_data), and the out_format.  It will try to overwrite anything
# that exists.  Returns 1 for success, undef for an error, and 0 for a
# close error (or if you give some bad data) (prints a message for any error).
sub save_data($$$){
	my($name,$channels,$out_format)=@_;
	my($OUT);

	# NOTE: to avoid confusion with the close and open methods above, close
	# and open will explicitly mention the CORE class.

	# Initialize stuff.
	# Try to get a date, somehow.
	my $datestrg;
	defined($datestrg=`date`) || ($datestrg=localtime);chomp $datestrg;
	# Make an array of what channels to include.
	# Just to be paranoid, make sure all the required elements exits.
	my @chan=();
	foreach(sort keys(%{$channels})){
		defined($channels->{$_}) || next;
		if(!(defined($channels->{$_}->{"tref"}) &&
		     defined($channels->{$_}->{"tinc"}) &&
		     defined($channels->{$_}->{"torig"}) &&
		     defined($channels->{$_}->{"yref"}) &&
		     defined($channels->{$_}->{"yinc"}) &&
		     defined($channels->{$_}->{"yorig"}) &&
		     defined($channels->{$_}->{"data"}))){
			print "WARNING: Skipping $_ because it seems to be malformed.\n";
			print "         This is an error with this script, and not the scope.\n";
			next;
		}
		push @chan,$_;
	}
	if($#chan<0){
		print "Found no valid channels.  Possible causes:\n";
		print "  1) You did not enable any channels for reading.\n";
		print "  2) This script has a bug.\n";
		print "Will not write data to $name.\n";
		return 0;
	}

	# Open file.
	if(!CORE::open($OUT,">",$name)){
		print "Could not open $name for writing.\n";return undef;
	}

	# Prepare and write the header.
	my $header="";
	$header.="# Oscilloscope data\n# Created: $datestrg\n";
	$header.="# Measured value M, in terms of scope units m, is: (m-ref)*inc+orig.\n";
	my $a;my $b;my $line="";for($b=0;$b<=$#chan;++$b){
		$header.="# $chan[$b] (time): ref=".$channels->{$chan[$b]}->{"tref"}.	
		                            " inc=".$channels->{$chan[$b]}->{"tinc"}.
		                           " orig=".$channels->{$chan[$b]}->{"torig"}."\n";
		$header.="# $chan[$b] (data): ref=".$channels->{$chan[$b]}->{"yref"}.	
		                            " inc=".$channels->{$chan[$b]}->{"yinc"}.
		                           " orig=".$channels->{$chan[$b]}->{"yorig"}."\n";
		if($b>0){$line.="\t";}
		if($out_format==0){$line.="Time ($chan[$b],scope units)\tData ($chan[$b],scope units)";}
		else{$line.="Time ($chan[$b],measured)\tData ($chan[$b],measured)";}
	}
	$header.="# ".$line."\n";
	if(!(print $OUT $header)){
		print "Error writing header to $name.\n";CORE::close($OUT);return undef;
	}

# This is how I used to write the data.  It was slower than what I do now:
#	# Assemble and write the data.
#	my $done;for($a=0;;++$a){
#		$line="";$done=1;for($b=0;$b<=$#chan;++$b){
#			if($b>0){$line.="\t";}
#			if($out_format==0){$line.=$a;}
#			else{$line.=($a-$channels->{$chan[$b]}{"tref"})*$channels->{$chan[$b]}{"tinc"}+$channels->{$chan[$b]}{"torig"};}
#			if(defined($channels->{$chan[$b]}{"data"}[$a])){
#				$done=0;$line.="\t";
#				if($out_format==0){$line.=$channels->{$chan[$b]}{"data"}[$a];}
#				else{$line.=($channels->{$chan[$b]}{"data"}[$a]-$channels->{$chan[$b]}{"tref"})*$channels->{$chan[$b]}{"tinc"}+$channels->{$chan[$b]}{"torig"};}
#			}else{$line.="\t";}
#		}
#		if($done){last;}
#		if(!(print $OUT "$line\n")){
#			print "Error writing data line $a to $name.\n";CORE::close($OUT);
#			return undef;
#		}
#	}
#	print "Wrote $a data points to $name.\n";

	# Assemble and write the data.
	# For speed, process the data in columns, and store in an array, where
	# each element stands for a line, and is a reference to an array for
	# every channel in that row.  This cuts down massively on the number of
	# references Perl has to trace through.  The speedup over doing things by
	# lines was moderate (factor of 2 or so).
	local $,="\t";
	my @data;for($b=0;$b<=$#chan;++$b){
		my $tref=$channels->{$chan[$b]}{"tref"};
		my $tinc=$channels->{$chan[$b]}{"tinc"};
		my $torig=$channels->{$chan[$b]}{"torig"};
		my $yref=$channels->{$chan[$b]}{"yref"};
		my $yinc=$channels->{$chan[$b]}{"yinc"};
		my $yorig=$channels->{$chan[$b]}{"yorig"};
		$a=0;foreach my $y (@{$channels->{$chan[$b]}{"data"}}){
			if(!defined($data[$a])){$data[$a]=[];}
			if($out_format==0){
				$data[$a][$b]="$a$,$y";
			}else{
				$data[$a][$b]=(($a-$tref)*$tinc+$torig) . $, .  (($y-$yref)*$yinc+$yorig);
			}
			++$a;
		}
	}
	$a=0;foreach(@data){
		if(!((print $OUT @{$_}) && (print $OUT "\n"))){
			print "Error writing data line $a to $name.\n";CORE::close($OUT);
			return undef;
		}
		++$a;
	}
	print "Wrote $a data points to $name.\n";

	# Close the file.
	if(!CORE::close($OUT)){
		print "Error closing $name.\n";return 0;
	}

	# Success:
	return 1;
}


# End of Agilent_Scope class.
}

return 1;
