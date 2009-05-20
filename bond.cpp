#include "bond.h"

int Bond::num_bonds = 0;

Bond::Bond () : bond(unbonded) {

return;
}

Bond::Bond (double length, bondtype btype) : bondlength(length), bond(btype) {

return;
}

Bond::~Bond () {

	--Bond::num_bonds;

return;
}

void Bond::SetBondType (const bondtype btype) {

	bond = btype;

	if (bond != unbonded)
		++Bond::num_bonds;

return;
}
