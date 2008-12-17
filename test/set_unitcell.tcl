# vmd tcl procedure: set the unitcell parameters for the whole trajectory
#
# Time-stamp: <akohlmey 02.07.2003 10:38:17 yello.theochem.ruhr-uni-bochum.de>
#
# Copyright (c) 2003 by <Axel.Kohlmeyer@theochem.ruhr-uni-bochum.de>
#

# arguments: 
#  molid = molecule id where the unitcell is added to (default top)

proc set_unitcell {a b c {molid top} {alpha 90.0} {beta 90.0} {gamma 90.0}} {

    if {![string compare $molid top]} {
        set molid [molinfo top]
    }

    set n   [molinfo $molid get numframes]

    for {set i 0} {$i < $n} {incr i} {
        molinfo $molid set frame $i
        molinfo $molid set {a b c alpha beta gamma} \
            [list $a $b $c $alpha $beta $gamma]
    }
}

############################################################
# Local Variables:
# mode: tcl
# time-stamp-format: "%u %02d.%02m.%y %02H:%02M:%02S %s"
# End:
############################################################
