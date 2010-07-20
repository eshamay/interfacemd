#include "carbonchain.h"

int CarbonChain::numCarbonChains = 0;

CarbonChain::CarbonChain (const int numCarbons) 
  : Molecule (), _carbons(Atom_ptr_vec(numCarbons, (AtomPtr)NULL))
{
  ++numCarbonChains;
}

CarbonChain::~CarbonChain () {
  numCarbonChains--;
}

CarbonChain::CarbonChain (const Molecule& molecule) 
  : Molecule(molecule)
{
  ++numCarbonChains;
}

// Returns a vector that points from the molecule's center of mass to the C10 atom - gives a rough approximation of the molecule's long-axis
VecR CarbonChain::Vector_CoM_To_End () {
	// First we set the atoms, update the center of mass, etc.
	this->SetAtoms ();
	this->UpdateCenterOfMass ();

	// the last carbon in the chain
	Atom * LastCarbon = _carbons[_carbons.size() - 1];
	// the axis between the center of mass and the last chain-carbon
	VecR axis (_centerofmass - LastCarbon->Position());
	return (axis);
}


// Runs through the list of atoms in the system and returns all the carbons
Atom_ptr_vec Carbons() {

  Atom_ptr_vec carbons;

  for (Atom_it ai = this->begin(); ai != this->end(); ai++)
  {
    // check if the atom has a name starting with 'C', and then followed by numbers
    if ((*ai)->Element() == Atom::C)
      carbons.push_back (*ai);
  }

  return carbons;
}
