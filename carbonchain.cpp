#include "carbonchain.h"

int CarbonChain::numCarbonChains = 0;

CarbonChain::CarbonChain (int numCarbons) 
  : Molecule (), _carbons(vector<Atom *>(numCarbons, (Atom *)NULL))
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

/*
// sets the _carbons by name for easier access
void CarbonChain::SetCarbons () {

  for (Atom_it ai = this->begin(); ai != this->end(); ai++)
  {
    std::string name = (*ai)->Name();
    if (name[0] == "C" && std::isdigit(name[1]))
      _carbons[FindDigits(name)] = *ai;
  }

  return;
}
*/
