#ifndef DECANE_H_
#define DECANE_H_

#include "carbonchain.h"
#include "utility.h"

// A water class to add a few functions for dealing with water molecules specifically.
class Decane: public CarbonChain {

	public:
		Decane ();	// a default constructor
		~Decane ();	// a destructor
		Decane (const Molecule& molecule);	// copy constructor for casting from a molecule

		static int numDecanes;			// total number of waters in the system

};

typedef std::vector<Decane *> Decane_ptr_vec;
typedef std::vector<Decane> Decane_vec;

#endif
