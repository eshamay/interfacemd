#!/usr/bin/perl -w 

# Author: JJT
# This provides my low-level TCP socket interface that can read and write
# through sockets with built-in timeouts.

{ package Low_Socket;
use strict;
use warnings;
use Socket;


# Constructor:
#   open -- the constructor.  Give it the host and port.
#           Call it directly (Low_Socket::open(), not Low_Socket->open()).

# Methods:
#   read -- reads.  Give it a timeout (undef for infinite).
#   write -- writes.  Give it the stuff to write.
#   close -- closes the socket.
#   DESTROY -- the destructor.  Just calls close above.

# The data used by this module is stored in a hash referenced by the class
# instance.  The keys are:
#   CONN -- the socket handle.
# Every instance of this class is assumed to be an open connection.  The
# connections are closed when the destructor is called.  If a connection is
# closed from the other end, the read and write methods should still return
# in a timely manner, but may give errors.  You can close a connection by
# calling the destructor (directly or by running "undef"), or by calling
# close.  I do not know of a way to check whether a socket is open, save to
# try writing something and seeing if that many bytes actually got written.
# You cannot reopen a socket -- you have to create a new one.
# NOTE:  The methods check for some internal consistency of the internal
# data, but it is possible to cause errors by changing the internal data
# from outside those methods, so do that at your own risk.


# Opens the connection, given the host and port.
# Returns an object for success, undef for error.
sub open($$){
	my($host,$port)=@_;
	my(%self,$iaddr,$paddr,$proto);

	# Make sure this was not called like Low_Socket->open():
	ref($host) && return undef;

	# Some checks.
	if(!$host){return undef;}
	if(!$port){return undef;}
	if($port=~/\D/){$port=getservbyname($port,"tcp") || return undef;}

	# Get some important stuff.
	($proto=getprotobyname("tcp")) || return undef;
	($iaddr=inet_aton($host)) || return undef;
	($paddr=sockaddr_in($port,$iaddr)) || return undef;

	# Make the socket and connection.
	socket($self{"CONN"},PF_INET,SOCK_STREAM,$proto) || return undef;
	connect($self{"CONN"},$paddr) || return undef;

	# Done.  Make it binary and autoflushing, and then return.
	my $oldFH=select($self{"CONN"});$|=1;select($oldFH);
	binmode($self{"CONN"});bless(\%self);
}


# Sends a string over the connection.
# Don't buffer anything, so use syswrite.
# Returns 0 for success, undef for a write error.
sub write($){
	my($this,$strg)=@_;

	# This needs to be a connected instance.
	if(!defined(ref($this))){return undef;}

	my($s,$size);
	if(!$strg){return 0;}
	# Check for an open connection, and whether CONN seems to be open.
	if(!defined(fileno($this->{"CONN"}))){return undef;}
	$size=0;while($s=syswrite($this->{"CONN"},substr($strg,$size))){
		# This may be superfluous:
		(defined($s)) || last;
		$size+=$s;($size>=length($strg)) && last;
	}
	($size<length($strg)) && return undef;
	return 0;
}


# Reads data from the connection, unbuffered.
# Returns what was read, which will be an empty string if nothing happened,
# and undef if there was an error.  A likely error is that the connection
# was closed from the other end.  The only errors are if the socket is not
# open, this was called incorrectly, or the connection was closed.
# Send a time argument for how long to wait (undef for infinite).
sub read($){
	my($this,$wait)=@_;

	# This needs to be a connected instance.
	if(!defined(ref($this))){return undef;}

	my $bits='';
	my $strg='';
	# Check for an open connection, and whether CONN seems to be open.
	if(!defined(fileno($this->{"CONN"}))){return undef;}
	vec($bits,fileno($this->{"CONN"}),1)=1;
	select($bits,undef,undef,$wait);
	if(vec($bits,fileno($this->{"CONN"}),1)){
		sysread($this->{"CONN"},$strg,65536);
		# If select says there was something to read, but there wasn't, that
		# seems to be an indication that the connection was closed.
		if(length($strg)==0){return undef;}
	}
	return $strg;
}


# Closes the connection.
# Doesn't return anything.
sub close(){
	my($this)=@_;
	if(defined(ref($this)) && defined($this->{"CONN"})){
		CORE::close($this->{"CONN"});
	}
}


# The destructor.  Just calls close.
sub DESTROY(){
	# Use the '&' to avoid confusion with the CORE::close() call.
	&close(@_);
}


# End of Low_Socket class.
}

return 1;
