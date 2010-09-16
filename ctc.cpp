#include "ctc.h"

CarbonTetrachloride::CarbonTetrachloride () {
	this->Rename("ctc");
	_moltype = Molecule::CTC;
}


CarbonTetrachloride::CarbonTetrachloride (const Molecule& mol)
	:
		Molecule (mol)
{
	this->Rename("ctc");
	_moltype = Molecule::CTC;
}


CarbonTetrachloride::CarbonTetrachloride (const MolPtr& mol)
	:
		Molecule (*mol)
{
	this->Rename("ctc");
	_moltype = Molecule::CTC;
}

void CarbonTetrachloride::SetAtoms () {

	/*
	_s = this->GetAtom(Atom::S);
	_o1 = (AtomPtr)NULL; _o2 = (AtomPtr)NULL;

	for (Atom_it it = this->begin(); it != this->end(); it++) {
		if ((*it)->Element() == Atom::O) {
			if (_o1 == NULL)
				_o1 = *it;
			else
				_o2 = *it;
		}
	}

	this->_so1 = this->_o1->Position() - this->_s->Position();
	this->_ctc = this->_o2->Position() - this->_s->Position();
	*/
	return;
}
