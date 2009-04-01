#!/usr/bin/perl -w


# An add-on to the Experiment library, giving Soloist support.
# Requires the Experiment library.
# It also exports a few function calls.


BEGIN {
package Experiment_Soloist;

use strict;
use warnings;
use Exporter qw(import);
use Experiment ();


# GLOBALS #############################


# The port to connect to on a Soloist:
my $soloist_port=24;
# The pulse time for the Soloist triggers (milliseconds):
my $soloist_trigger=5;


# This is a list of the known Soloist boxes.  Each one has a name (used to
# connect to it) and an output on a digital box.  That output on the
# digital box is used to trigger motion on the Soloist.
# This is actually a hash.  The keys are the names of the Soloists, and the
# values are hash refs.  Those hashes have two keys: box (the digital box
# name) and out (the number of the output on the box).
# A third key is optional: program, if it exists, is the last program sent
# to the box.
my %soloist=();


# INTERNAL ROUTINES ###################


# Given a list of actions (described below -- these had better be sorted by
# time), this forms a hash.  The keys are all the soloist boxes referenced
# (by name), and the values are hash references, of which the important
# keys are 'data', which contain the string with which the box should be
# programmed.  A more complete description of the hash is in the comments
# preceding the Experiment::program_boxes() routine (the box_data data
# type).
# Returns the hash reference if successful, and undef for error.  
# You can go ahead and give all the actions here -- it only pays attention
# to the Soloist_Action s.
# Warns about various likely-to-be-problems, but can't catch everything.
# It does not know for sure how long a move on the Soloist will take, so it
# cannot catch when you are interrupting a move.  It will catch when the
# triggers for one particular Soloist overlap (and will return undef), but
# does not check if triggers for different Soloists, which may run off the
# same trigger, overlap.
sub soloist_box_data(@){
	my @action_list=@_;
	my($action,$box,$delay);
	my %soloist_box_data=();
	foreach $action (@action_list){
		# Ignore things not meant for this routine.
		if(ref($action) ne "Soloist_Action"){next;}
		$box=$action->{'box'};
		if(!defined($soloist_box_data{$box})){
			# Initialize data to "\n", which will empty the Soloist buffer of
			# commands.
			$soloist_box_data{$box}={
				'action_type'=>'Soloist_Action',
				'data'=>"\n",
				'time'=>$action->{'time'},
			};
		}else{
			# Check if we're too soon after the last trigger.  I'm currently in a
			# paranoid mood, where I use the negation of '>' instead of '<=' so
			# that this test will pass and a warning will happen if something is
			# not numeric.
			if(!($action->{'time'}>($soloist_box_data{$box}{'time'}+$soloist_trigger))){
				carp("soloist_action(): trigger for $box overlaps previous at time $action->{'time'}.\n");
				return undef;
			}
			$soloist_box_data{$box}{'time'}=$action->{'time'};
		}
		# The triggers should already be in the event queue, so don't bother
		# with that.  Just write the program.  Of course, the program is
		# already written, really.  I just need to concatenate the lines.
		# Two cases -- if the position is undefined, give a -1 to disable the
		# air cart.  Otherwise, give the position and speed.
		foreach my $ref (@{$action->{'points'}}){
			if(!defined($ref->[0])){
				$soloist_box_data{$box}{'data'}.="-1\n";
			}else{
				$soloist_box_data{$box}{'data'}.="$ref->[0] $ref->[1]\n";
			}
		}
	}
	return \%soloist_box_data;
}


# Given a box_data structure (one particular box, described before the
# program_boxes() routine), this interprets the data that was sent, and
# responds to any requests.  The second argument is the most recently read
# stuff from the box.  The third argument is the box name, for more useful
# error messages (this prints error messages).  Returns undef for error,
# negative if box said it was full, 0 if things are still going, and
# postive if we have finished sending data and the connection may be
# closed.
# This particular version is specific to the Soloist.
sub respond_to_soloist($$$){
	my($box_data,$data,$box)=@_;
	# Append data to req:
	if(!defined($box_data->{'req'})){$box_data->{'req'}=$data;}
	else{$box_data->{'req'}.=$data;}
	# Make sure box type is correct:
	if($box_data->{'action_type'} ne 'Soloist_Action'){
		carp("respond_to_soloist(): unknown action type ".$box_data->{'action_type'}." for $box.\n");
		return undef;
	}
	# Send data only if this is the first time, or if there is something to
	# respond to (Soloist's responses end in '\n'):
	if(length($box_data->{'sent'})==0 || $box_data->{'req'}=~/\n$/o){
		# If the Soloist has said anything, that should be a number telling us
		# how many lines left it has in its buffer.  If we are here, it's
		# because there is still something left to send, so a 0 means we are in
		# trouble.
		if(length($box_data->{'req'})>0){
			my $left=$box_data->{'req'};
			$left=~s/[\r\n]+$//o;
			if(!($left=~/^\d+$/o)){
				carp("respond_to_soloist(): bad request from $box.\n");
				next;
			}
			# Return if full.
			if($left==0){return -1;}
		}
		# Send one line of the program.
		my($send,$remainder)=($box_data->{'data'}=~/^(.*?)[\r\n]+([^\r\n].*)$/os);
		$send="$send\n";$box_data->{'data'}=$remainder;
		if(!defined($send) || $send eq ''){
			carp("respond_to_soloist(): data seems corrupted for $box.\n");
			return undef;
		}
		# For writing to a network, writing should be limited by network
		# speed and stuff, so there's no reason not to let this block.
		# Make a quick check that the socket seems open.
		if(!defined(fileno($box_data->{'sock'}))){
			carp("respond_to_soloist(): socket for $box does not appear to be open.\n");
			return undef;
		}
		while(my $s=syswrite($box_data->{'sock'},substr($box_data->{'data'},0,$total))){
			if(!defined($s)){
				carp("respond_to_soloist(): problem sending data to $box ($!).\n");
				return undef;
			}
			# Move the written portion from send to sent.
			$box_data->{'sent'}.=substr($send,0,$s);
			$send=substr($send,$s);
			if(length($send)==0){last;}
		}
		# We have responded to this response:
		$box_data->{'req'}='';
		# We are done here -- signal if we are done completely
		if(length($box_data->{'data'})==0){return 1;}
	}
	# We are done for now, but not done sending data.
	return 0;
}


# ACTION/EVENT MANIPULATION ###########


# Use this to declare Soloist boxes.  Basically, give the box name of the
# Soloist, the box name of the digital box that triggers it, and the
# channel number on that digital box controller.  These routines will test
# if you are attempting to trigger a Soloist too often, but have no test
# for whether you use the same channel on the same box to trigger multiple
# Soloists.  Using the same channel on the same box to trigger more than
# one thing, even if they should always be simultaneous, is frowned upon.
# It should be OK if the Soloists are always given the same times, but if
# you have one move when the other doesn't, that will screw things up.
# Don't give the wrong names (bad things happen if you give an analog box
# name for the trigger box, for instance).
# Returns undef for error, and true otherwise.
# Prints a warning if the given soloist name has been used, but returns
# true and overwrites the old one nevertheless.
sub declare_soloist($$$){
	my($box,$dig_box,$dig_out)=@_;
	if($dig_out<0 || $dig_out>31)){
		carp("declare_soloist(): channel $dig_out is out of range.\n");
		return undef;
	}
	if(defined($soloist{$box})){
		carp("declare_soloist(): $box has already been declared.  Will overwrite.\n");
	}
	$soloist{$box}={'box'=>$dig_box,'out'=>$dig_out};
	return 1;
}


# This returns something like a Soloist_Action.  Actually, it returns an
# Event, which contains Digital_Action's to trigger the Soloist, and a
# Soloist_Action.  Every Soloist_Action is just a list of [position,speed]
# pairs.  These will get sent to the Soloist directly.  Between actions,
# there will be a command to tell the Soloist to wait for a digital
# trigger.
# WARNING:  It is up to the user to make sure that the Soloist can complete
# the motions before the next trigger.  There are no checks for that here.
# The Soloist_Action type is fairly simple:
# The keys are:
#   'time'   => the time this action takes place (in milliseconds)
#   'box'    => the name of the soloist box to program
#   'points' => a list reference of points.  Each point is a ref to a
#               2-element list.  The first is the position, the second
#               is the speed.  The Soloist will attempt to run through all
#               these in the order given.  Subsequent actions will require
#               a digital trigger (included with the action at not cost).

# This routine takes the following:
# A box name.  This is the name of the Soloist controller.  It must be a
# declared box (see declare_soloist()).
# A bunch (possibly 0, but why?) of [position,speed] pairs of where to
# move, and how fast.  Positions are in mm, and speeds are in mm/s.
# Make the position undef to disable the air cart.  The speed will be
# ignored in that case.
# Returns an Event, containing the Digital_Actions for the trigger and the
# Soloist_Action.
# Will warn you if you have bad arguments.
sub soloist_action($@){
	my($box)=shift(@_);
	# Initialization.
	# All actions are defined to be at time 0 -- add them to an event to set
	# the time.
	if(!defined($soloist{$box})){
		carp("soloist_action(): $box is not a declared Soloist box.\n");
		return undef;
	}
	my $sol_ref={'time'=>0,'box'=>$box,'points'=>[]};
	my $sol_dig=$soloist{$box}{'box'};
	my $sol_out=$soloist{$box}{'out'};
	bless($sol_ref,"Soloist_Action");
	# Add the waypoints, checking as you go.
	foreach $pt (@_){
		if(ref($pt) ne "ARRAY" || scalar(@$pt)!=2){
			carp("soloist_action(): needs references to 2-element arrays only.\n");
			return undef;
		}
		# Sanity checks on the given points.  There is range-checking here,
		# too.
		if(!defined($pt->[0])){
			# Make both undef.  That means disable the air cart.
			push(@{$sol_ref->{'points'}},[undef,undef]);
			next;
		}
		if(!defined($pt->[1])){
			carp("soloist_action(): 2-element array not entirely defined.\n");
			return undef;
		}
		if($pt->[0]<0 || $pt->[0]>100)){
			carp("soloist_action(): position $pt->[0] is not within range.\n");
			return undef;
		}
		if($pt->[1]<0.1 || $pt->[1]>50)){
			carp("soloist_action(): speed $pt->[1] is not considered reasonable.\n");
			return undef;
		}
		# Make a copy of the array, so that the internals here won't change.
		push(@{$sol_ref->{'points'}},[$pt[0],$pt[1]]);
	}
	# Return the full event.
	return event(0,$sol_ref,digital_action($sol_dig,[$sol_out],[]),
	             event($soloist_trigger,digital_action($sol_dig,[],[$sol_out])));
}


# FIXME -- May need to dump the stream handle.

# FIXME -- Should change Experiment library so it can handle add-ons.  A
# hash, where keys are action types (Digital_Action, Analog_Action,
# Soloist_Action), pointing to hash refs.  In each, there should be
# pointers to the subroutines to call.  The box_data routine, the
# respond_to_box routine, something to copy actions (for the event
# routine), etc.  program_boxes() and event() will need to be changed.
# Others, too, perhaps.
# Should change the "type" argument in the box_data structure to be the
# action name (type -> action_type).  That will also require changing
# respond_to_box().
# Make other changes, like max number of sockets, to program_boxes, and
# change some error messages (box $box is probably redundant).  Have it
# call the respond_to_box-type routine on opening the socket with null
# response -- modify respond_to_box() to do nothing in that case.  The
# Soloist version needs it.

# FIXME -- how to handle pinging?


# EXPORT STUFF ########################


# FIXME -- should export some things, methinks.

# FIXME -- should actually add self to Experiment::%Action.

# If you export only subroutines, it is faster to leave off the '&' on each
# one (from the Exporter documentation).  If no type is given, the '&' is
# assumed.
our @EXPORT=qw();
our @EXPORT_OK=qw();


}
1;
