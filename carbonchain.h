#ifndef CARBONCHAIN_H_
#define CARBONCHAIN_H_

#include "molecule.h"
#include "utility.h"

class CarbonChain : public molecule {
	protected:
		vector<Atom *> _carbons;		/* An ordered listing of all the carbons in
										   the molecule */
	public:
		CarbonChain ();	// a default constructor
		virtual ~CarbonChain ();	// a destructor
		CarbonChain (const Molecule& molecule);	// copy constructor for casting from a molecule

		static int numCarbonChains;			// total number of waters in the system

		// Functions for analysis
		virtual void SetAtoms ();

		vector<Atom *>& Carbons () { return (_carbons); }
		Atom * Carbon (int index) { return (_carbons[index]); }

		VecR Vector_CoM_To_End ();
};

typedef std::vector<CarbonChain *> CarbonChain_ptr_vec;
typedef std::vector<CarbonChain> CarbonChain_vec;
