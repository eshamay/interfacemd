#include "oh.h"


	int Hydroxide::numHydroxides = 0;

	Hydroxide::Hydroxide () : Molecule()
	{
		this->Rename("oh");
		this->_moltype = Molecule::OH;
		++numHydroxides;
	}

	Hydroxide::~Hydroxide () {
		--numHydroxides;
	}

	void Hydroxide::SetAtoms () {
		this->FixAtoms();
		_o = this->GetAtom("O");
		_h = this->GetAtom("H");

		_oh = _h->Position() - _o->Position();

		return;
	}

