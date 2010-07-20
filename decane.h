#ifndef DECANE_H_
#define DECANE_H_

#include "alkane.h"

// A water class to add a few functions for dealing with water molecules specifically.
class Decane: public Alkane {

	public:
		Decane ();	// a default constructor
		~Decane ();	// a destructor
		Decane (const Molecule& molecule);	// copy constructor for casting from a molecule

		static int numDecanes;			// total number of waters in the system

		void SetAtoms ();

};

typedef std::vector<Decane *> Decane_ptr_vec;
typedef std::vector<Decane> Decane_vec;

#endif
