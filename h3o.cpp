#include "h3o.h"

int Hydronium::numHydroniums = 0;

Hydronium::Hydronium () : Molecule()
{
	_name = "h3o";
	++numHydroniums;
}

Hydronium::~Hydronium () {
	--numHydroniums;
}

void Hydronium::SetAtoms () {

		// here's the hydrogen and nitrogen atoms
		_o = (*this)["O"];

		// now set the 3 hydrogens
		std::vector<Atom *> hydrogens;
		RUN (_atoms) {
			if (_atoms[i]->Name() != "H") continue;
			hydrogens.push_back (_atoms[i]);
		}

		_h1 = hydrogens[0];
		_h2 = hydrogens[1];
		_h3 = hydrogens[2];

		// while we're here we may as well also find the N-O bond vectors
		_oh1 = _h1->Position() - _o->Position();
		_oh2 = _h2->Position() - _o->Position();
		_oh3 = _h3->Position() - _o->Position();

		_set = true;

return;
}
