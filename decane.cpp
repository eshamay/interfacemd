#include "decane.h"

int Decane::numDecanes = 0;

Decane::Decane () 
  : Alkane (10)
{
  this->Rename("dec");
  _moltype = Molecule::DECANE;
  ++numDecanes;
}

Decane::~Decane () {
  numDecanes--;
}

Decane::Decane (const Molecule& molecule) 
  : Alkane(molecule)
{
  ++numDecanes;
}

void Decane::SetAtoms () {

  // Runs through the list of atoms in the system and returns all the carbons
  typedef const Atom::Element_t& (Atom::*fn)() const;
  md_utility::copy_if (this->begin(), this->end(), std::back_inserter(_carbons), md_utility::mem_fun_eq(&Atom::Element, Atom::C));

  /*
  // Go through and assign each carbon to the slot in the ordered list
  _carbons[0] = this->GetAtom("C1");
  _carbons[1] = this->GetAtom("C2");
  _carbons[2] = this->GetAtom("C3");
  _carbons[3] = this->GetAtom("C4");
  _carbons[4] = this->GetAtom("C5");
  _carbons[5] = this->GetAtom("C6");
  _carbons[6] = this->GetAtom("C7");
  _carbons[7] = this->GetAtom("C8");
  _carbons[8] = this->GetAtom("C9");
  _carbons[9] = this->GetAtom("C10");
  */
	
  return;
}
