#include "bond.h"

int Bond::num_bonds = 0;

Bond::Bond () :
	bond(unbonded)
{
	++Bond::num_bonds;
}

Bond::Bond (double length, bondtype btype) :
	bondlength(length),
	bond(btype)
{
	++Bond::num_bonds;
}

Bond::Bond (const Bond& oldBond) :
	bondlength(oldBond.bondlength),
	bond(oldBond.bond)
{
	++Bond::num_bonds;
}

Bond::~Bond ()
{
	--Bond::num_bonds;
}

void Bond::SetBondType (const bondtype btype) {

	bond = btype;

	if (bond != unbonded)
		++Bond::num_bonds;

return;
}
