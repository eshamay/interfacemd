#ifndef CARBONCHAIN_H_
#define CARBONCHAIN_H_

#include "molecule.h"

class CarbonChain : public Molecule {

  protected:
    vector<Atom *> _carbons;			/* An ordered listing of all the carbons in
						   the molecule */
  public:
    CarbonChain (int numCarbons);			// a default constructor
    virtual ~CarbonChain ();
    CarbonChain (const Molecule& molecule);		// copy constructor for casting from a molecule

    static int numCarbonChains;			// total number of carbon chains in the system

    // Functions for analysis
    virtual void SetAtoms () = 0;
    void SetCarbons ();

    Atom_ptr_vec& Carbons () { return (_carbons); }
    Atom * Carbon (int index) { return (_carbons[index]); }

    VecR Vector_CoM_To_End ();
};

typedef std::vector<CarbonChain *> CarbonChain_ptr_vec;
typedef std::vector<CarbonChain> CarbonChain_vec;

#endif
