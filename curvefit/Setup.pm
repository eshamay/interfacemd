#!/usr/bin/perl -w


# A library containing stuff pertinent to the physical setup.

BEGIN {
package Setup;

use strict;
use warnings;
use Experiment;
use Exporter;
# This includes the import method of Exporter as part of this package
# (actually, includes all of Exporter as part of this), so that it gets
# called when you use this package (and exports the stuff listed in
# @EXPORT and so on).
our @ISA=qw(Exporter);


# Standard values #####################
# Here is where to store the values get_mot.pl should set things.  Values
# used in other scripts should perhaps use these as baseline values.

# Values for the repump amplitude and frequency.
# The repump intensity has been changed to full computer control, so the
# knob should be at 0.  Use 6 to get a decent MOT.  That's not optimized --
# the MOT still looks good at 3.  1.5 seems to be a good starting point for
# searching for dipole trap loading.
# For the repump intensity, the knob was around 11:00.
our $REPUMP_FREQ=0;
our $REPUMP_INTENSITY=6;

# These are the values for the MOT frequency and MOT2 intensity:
our $MOT_FREQ=7.25;
our $MOT2_INTENSITY=7.5;
# POSSIBLY TEMPORARY -- We have a different AOM in the MOT2 setup, and it
# takes a different intensity than we used to use.  It was swapped in on
# 2007-06-05 (an Isomet -- see lab notebook):
#$MOT2_INTENSITY=5.5;
# We went back to the IntraAction AOM on 2007-07-17, so I just commented
# the above line out.  The previous value was still a good maximum.

# These are the coil values for the MOT2 Helmholtz coils:
our $MOT2_COIL_X=-0.14;
our $MOT2_COIL_Y=-0.15;
our $MOT2_COIL_Z=-0.36;
# We swapped some fibers on the MOT setup, and so now have different coil
# settings:
$MOT2_COIL_X=0.02;
$MOT2_COIL_Y=-0.44;
$MOT2_COIL_Z=-0.42;

# These are offsets for the MOT2 AH coils (master and slave):
our $MOT2_AH_M=0;
our $MOT2_AH_S=0;


# Standard actions ####################


# Camera pulse, and reset it.
# Old version, before relay got weak:
#our $cam_trg        =digital_action('digital1',[16,0],[]);
#our $cam_rst        =digital_action('digital1',[],[16,0]);
# New version:
our $cam_trg        =analog_action_set('analog1',0,5);
our $cam_rst        =analog_action_set('analog1',0,0);
# Repump beam on and off.
our $repump_beam_on =digital_action('digital1',[],[17,1]);
our $repump_beam_off=digital_action('digital1',[17,1],[]);
# MOT1 beam on and off.
# IMPORTANT:  This triggers the MOT1 shutter, too.  However, the timing of
# that has not been tested, so don't expect this to turn on and off exactly
# when you say (the AOMs will turn it off right away, and the shutter will
# be a bit delayed, but the shutter won't let the beam turn back on right
# away).
our $mot1_beam_on   =digital_action('digital1',[26,10],[18,2]);
our $mot1_beam_off  =digital_action('digital1',[18,2],[26,10]);
# MOT2 beam on and off.
our $mot2_beam_on   =digital_action('digital1',[],[19,3]);
our $mot2_beam_off  =digital_action('digital1',[19,3],[]);
# MOT1 coils on and off.
our $mot1_coil_on   =digital_action('digital1',[20,4],[]);
our $mot1_coil_off  =digital_action('digital1',[],[20,4]);
# MOT2 coils on and off.
our $mot2_coil_on   =digital_action('digital1',[21,5],[]);
our $mot2_coil_off  =digital_action('digital1',[],[21,5]);

# A trigger pulse to start the other boards from analog1:
our $trigger_boards =event(0,
                           event(0,analog_action_set('analog1',15,5)),
                           event(100,analog_action_set('analog1',15,0)));

# These are probably temporary, and this pin will probably be used later.
# For the one-way barrier, we have two beams -- the barrier beam itself
# (tuned to an Rb 85 line, or somewhere between the repump and MOT
# transitions of Rb 87), which is controlled by a mechanical shutter (don't
# need precise timing for that), and the repump part of the barrier, which
# is the 0 order of the repump beam.  That means it is off by 80Mhz, so we
# use another 80MHz AOM to turn that on and off.  This turns that AOM on
# and off.  When the oneway AOM is on, the oneway beam will be on.  When
# the repump AOM is off, there is more in the 0 order, so the oneway beam
# will be brighter:
our $oneway_repump_on =digital_action('digital1',[],[22,6]);
our $oneway_repump_off=digital_action('digital1',[22,6],[]);
# And these trigger the shutter for the main beam.  oneway_beam_on turns
# the beam on (turns the shutter off):
our $oneway_beam_on =digital_action('digital1',[],[28,12]);
our $oneway_beam_off=digital_action('digital1',[28,12],[]);

# These turn the dipole beam on and off via shutter:
our $dipole_beam_on =digital_action('digital1',[],[23,7]);
our $dipole_beam_off=digital_action('digital1',[23,7],[]);


# Routines that give you an action ####


# Set the repump intensity to the given value:
# TEMPORARAY?  Channel 14 doesn't seem to work, so moved to 2.
sub set_repump_intensity($){
	if($_[0]<0){die("set_repump_intensity():  Will not set negative.\n");}
	return analog_action_set('analog1',2,$_[0]);
}
# Ramp the repump intensity from the last set to the given value:
sub ramp_repump_intensity($){
	if($_[0]<0){die("ramp_repump_intensity():  Will not set negative.\n");}
	return analog_action_ramp('analog1',2,$_[0]);
}

# Set the MOT master detuning to the given value:
sub set_mot_freq($){
	if($_[0]<0){die("set_mot_freq():  Will not set negative.\n");}
	return analog_action_set('analog1',13,$_[0]);
}
# Ramp the MOT master detuning from the last set to the given value:
sub ramp_mot_freq($){
	if($_[0]<0){die("ramp_mot_freq():  Will not set negative.\n");}
	return analog_action_ramp('analog1',13,$_[0]);
}

# Set the MOT1 intensity to the given value:
sub set_mot1_intensity($){
	if($_[0]<0){die("set_mot1_intensity():  Will not set negative.\n");}
	return analog_action_set('analog1',12,$_[0]);
}
# Ramp the MOT1 intensity from the last set to the given value:
sub ramp_mot1_intensity($){
	if($_[0]<0){die("ramp_mot1_intensity():  Will not set negative.\n");}
	return analog_action_ramp('analog1',12,$_[0]);
}

# Set the MOT2 intensity to the given value:
sub set_mot2_intensity($){
	if($_[0]<0){die("set_mot2_intensity():  Will not set negative.\n");}
	return analog_action_set('analog1',11,$_[0]);
}
# Ramp the MOT2 intensity from the last set to the given value:
sub ramp_mot2_intensity($){
	if($_[0]<0){die("ramp_mot2_intensity():  Will not set negative.\n");}
	return analog_action_ramp('analog1',11,$_[0]);
}

# Set/ramp the MOT2 AH master coil channel voltage:
sub set_mot2_AHM_coil($){
	return analog_action_set('analog1',10,$_[0]);
}
sub ramp_mot2_AHM_coil($){
	return analog_action_ramp('analog1',10,$_[0]);
}
# Set/ramp the MOT2 AH slave coil channel voltage:
sub set_mot2_AHS_coil($){
	return analog_action_set('analog1',9,$_[0]);
}
sub ramp_mot2_AHS_coil($){
	return analog_action_ramp('analog1',9,$_[0]);
}

# Set/ramp the MOT2 X coil channel voltage:
sub set_mot2_x_coil($){
	return analog_action_set('analog1',7,$_[0]);
}
sub ramp_mot2_x_coil($){
	return analog_action_ramp('analog1',7,$_[0]);
}
# Set/ramp the MOT2 Y coil channel voltage:
sub set_mot2_y_coil($){
	return analog_action_set('analog1',6,$_[0]);
}
sub ramp_mot2_y_coil($){
	return analog_action_ramp('analog1',6,$_[0]);
}
# Set/ramp the MOT2 Z coil channel voltage:
sub set_mot2_z_coil($){
	return analog_action_set('analog1',5,$_[0]);
}
sub ramp_mot2_z_coil($){
	return analog_action_ramp('analog1',5,$_[0]);
}


# TEMPORARY FIXES: ######################


# EXPORT STUFF ########################


our @EXPORT=qw($REPUMP_FREQ $REPUMP_INTENSITY $MOT_FREQ $MOT2_INTENSITY
               $MOT2_COIL_X $MOT2_COIL_Y $MOT2_COIL_Z
               $MOT2_AH_M $MOT2_AH_S
               $cam_trg $cam_rst $trigger_boards
               $repump_beam_on $repump_beam_off
               $mot1_beam_on $mot1_beam_off $mot2_beam_on $mot2_beam_off
               $mot2_coil_on $mot2_coil_off
               &set_repump_intensity &ramp_repump_intensity
               &set_mot_freq &ramp_mot_freq
               &set_mot1_intensity &ramp_mot1_intensity
               &set_mot2_intensity &ramp_mot2_intensity
               &set_mot2_AHM_coil &ramp_mot2_AHM_coil
               &set_mot2_AHS_coil &ramp_mot2_AHS_coil
               &set_mot2_x_coil &ramp_mot2_x_coil
               &set_mot2_y_coil &ramp_mot2_y_coil
               &set_mot2_z_coil &ramp_mot2_z_coil
               $dipole_beam_on $dipole_beam_off
               $oneway_beam_on $oneway_beam_off
               $oneway_repump_on $oneway_repump_off
              );

our @EXPORT_OK=();


}
1;
