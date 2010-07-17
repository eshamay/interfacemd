#include "h3o.h"

int Hydronium::numHydroniums = 0;

Hydronium::Hydronium () : Molecule()
{
	this->Rename("h3o");
	++numHydroniums;
}

Hydronium::~Hydronium () {
	--numHydroniums;
}

void Hydronium::SetAtoms () {

		// here's the hydrogen and nitrogen atoms
		_o = this->GetAtom(Atom::O);

		// now set the 3 hydrogens
		Atom_ptr_vec hydrogens;
		for (Atom_it it = this->begin(); it != this->end(); it++) {
			if (*it == _o) continue;
			hydrogens.push_back (*it);
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
