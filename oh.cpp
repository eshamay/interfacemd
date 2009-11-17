#include "oh.h"

int Hydroxide::numHydroxides = 0;

Hydroxide::Hydroxide () : Molecule()
{
	_name = "oh";
	++numHydroxides;
}

Hydroxide::~Hydroxide () {
	--numHydroxides;
}

void Hydroxide::SetAtoms () {
	_o = this->GetAtom("O");
	_h = this->GetAtom("H");

	_oh = _h->Position() - _o->Position();
	_set = true;

return;
}

