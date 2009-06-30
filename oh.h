#ifndef OH_H_
#define OH_H_

#include "h2o.h"

/******************************
 * Hydroxide (H3O+)
 * ****************************/
class Hydroxide: public Molecule {

protected:
	VecR _oh;				// The OH bond vector
	Atom *_o, *_h;			// pointers to the atoms for easy access

public:
	Hydroxide ();
	~Hydroxide ();
	Hydroxide (const Molecule& molecule);	// copy constructor for casting from a molecule

	static int numHydroxides;			// total number of waters in the system

	void SetAtoms ();					// set the _oh bond vector
	VecR const * OH () const { return &_oh; }
};

#endif
