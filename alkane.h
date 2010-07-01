#ifndef DECANE_H_
#define DECANE_H_

#include "carbonchain.h"

// A water class to add a few functions for dealing with water molecules specifically.
class Alkane: public CarbonChain {

  public:
    Alkane ();	// a default constructor
    ~Alkane ();	// a destructor
    Alkane (const Molecule& molecule);	// copy constructor for casting from a molecule

    static int numAlkanes;			// total number of waters in the system

    void SetAtoms ();

};

typedef std::vector<Alkane *> Alkane_ptr_vec;
typedef std::vector<Alkane> Alkane_vec;

#endif
