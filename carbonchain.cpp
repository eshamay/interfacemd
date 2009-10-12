#include "carbonchain.h"

int CarbonChain::numCarbonChains = 0;

CarbonChain::CarbonChain () 
  : Molecule () 
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
VecR CarbonChain::Vector_CoM_to_Atom (const Atom * atom) {
	// First we set the atoms, update the center of mass, etc.
	this->SetAtoms ();
	this->UpdateCenterOfMass ();

	VecR axis (_centerofmass.MinVector(atom));
	return (axis);
}
