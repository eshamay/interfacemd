#include "alkane.h"

int Alkane::numAlkanes = 0;

Alkane::Alkane () 
  : CarbonChain (10)
{
  _name = "dec";
  ++numAlkanes;
}

Alkane::~Alkane () {
  numAlkanes--;
}

Alkane::Alkane (const Molecule& molecule) 
  : CarbonChain(molecule)
{
  ++numAlkanes;
}

void Alkane::SetAtoms () {
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
	
  return;
}

// Runs through the list of atoms in the system and returns all the carbons
Atom_ptr_vec Carbons() {

  Atom_ptr_vec carbons;

  for (Atom_it ai = this->begin(); ai != this->end(); ai++)
  {
    // check if the atom has a name starting with 'C', and then followed by numbers
    name = (*ai)->Name();
    if (name[0] == "C" && std::isdigit(name[1]))
      carbons.push_back (*ai);
  }

  return carbons;
}

 
