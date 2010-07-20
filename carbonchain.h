#ifndef CARBONCHAIN_H_
#define CARBONCHAIN_H_

#include "molecule.h"

class CarbonChain : public Molecule {

  protected:
	Atom_ptr_vec _carbons;			/* An ordered listing of all the carbons in
						   the molecule */
  public:
    CarbonChain (const int numCarbons);			// a default constructor
    virtual ~CarbonChain ();
    CarbonChain (const Molecule& molecule);		// copy constructor for casting from a molecule

    static int numCarbonChains;			// total number of carbon chains in the system

    // Functions for analysis
    virtual void SetAtoms () = 0;
    void SetCarbons ();

	Atom_it carbons_begin () const { return _carbons.begin(); }
	Atom_it carbons_end () const { return _carbons.end(); }
    Atom_ptr_vec& Carbons () const { return (_carbons); }
    AtomPtr Carbon (const int index) const { return (_carbons[index]); }

    VecR Vector_CoM_To_End ();
};

typedef std::vector<CarbonChain *> CarbonChain_ptr_vec;
typedef std::vector<CarbonChain> CarbonChain_vec;

#endif
