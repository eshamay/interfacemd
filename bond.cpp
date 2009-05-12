#include "bond.h"

int Bond::num_bonds = 0;

Bond::Bond () : bond(unbonded) {

return;
}

Bond::Bond (double length) : bondlength(length) {

	SetBondType ();

return;
}

Bond::~Bond () {

	--Bond::num_bonds;

return;
}

void Bond::SetBondType () {

	bond = unbonded;

	// one type of bond is the O-H covalent
	if (bondlength <= OHBONDLENGTH)
		bond = ohbond;

	// Or an H-bond is formed!
	if (bondlength <= HBONDLENGTH and bondlength > OHBONDLENGTH)
		bond = hbond;

	if (bond != unbonded)
		++Bond::num_bonds;

return;
}
