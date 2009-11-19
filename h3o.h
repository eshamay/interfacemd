#ifndef H3O_H_
#define H3O_H_

#include "h2o.h"

/******************************
 * Hydronium (H3O+)
 * ****************************/
class Hydronium: public Molecule {

protected:
	VecR _oh1, _oh2, _oh3;				// Both of the OH vectors
	Atom *_o, *_h1, *_h2, *_h3;			// pointers to the atoms for easy access

public:
	Hydronium ();
	~Hydronium ();
	Hydronium (const Molecule& molecule);	// copy constructor for casting from a molecule

	static int numHydroniums;			// total number of waters in the system

	void SetAtoms ();
	VecR MolecularAxis () { return _z; }
	VecR const * OH1 () const { return &_oh1; }
	VecR const * OH2 () const { return &_oh2; }
	VecR const * OH3 () const { return &_oh3; }
};

#endif
