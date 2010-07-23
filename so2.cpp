#include "so2.h"

SulfurDioxide::SulfurDioxide () {
  this->Rename("so2");
  _moltype = Molecule::SO2;
}


SulfurDioxide::SulfurDioxide (const Molecule& mol)
  :
	Molecule (mol)
{
  this->Rename("so2");
  _moltype = Molecule::SO2;
}


SulfurDioxide::SulfurDioxide (const MolPtr& mol)
  :
	Molecule (*mol)
{
  this->Rename("so2");
  _moltype = Molecule::SO2;
}

void SulfurDioxide::SetAtoms () {

  _s = this->GetAtom("S");
  _o1 = (AtomPtr)NULL; _o2 = (AtomPtr)NULL;

  for (Atom_it it = this->begin(); it != this->end(); it++) {
	if ((*it)->Element() == Atom::O) {
	  if (_o1 == NULL)
		_o1 = *it;
	  else
		_o2 = *it;
	}
  }
  this->FixAtoms();
  return;
}
