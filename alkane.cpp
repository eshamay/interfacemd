#include "alkane.h"

int Alkane::numAlkanes = 0;

  Alkane::Alkane (const int numCarbons) 
: Molecule (), _carbons(Atom_ptr_vec(numCarbons, (AtomPtr)NULL))
{
  ++numAlkanes;
  this->Rename("alkane");
  _moltype = Molecule::ALKANE;
}

Alkane::~Alkane () {
  numAlkanes--;
}

  Alkane::Alkane (const Molecule& molecule) 
: Molecule(molecule)
{
  ++numAlkanes;
}

// Returns a vector that points from the molecule's center of mass to the C10 atom - gives a rough approximation of the molecule's long-axis
VecR Alkane::Vector_CoM_To_End () {
  // First we set the atoms, update the center of mass, etc.
  this->SetAtoms ();
  this->UpdateCenterOfMass ();

  // the last carbon in the chain
  Atom * LastCarbon = _carbons[_carbons.size() - 1];
  // the axis between the center of mass and the last chain-carbon
  VecR axis (_centerofmass - LastCarbon->Position());
  return (axis);
}


