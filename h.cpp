#include "h.h"

int Proton::numProtons = 0;

Proton::Proton () : Molecule()
{
	this->Rename("h+");
	_moltype = Molecule::H;
	++numProtons;
}

Proton::~Proton () {
	--numProtons;
}

void Proton::SetAtoms () {
	_h = this->GetAtom("H");
	this->FixAtoms();
return;
}

