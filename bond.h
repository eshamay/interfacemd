#ifndef BOND_H_
#define BOND_H_

#include <map>
#include <vector>
#include <math.h>
#include <string>

const double OHBONDLENGTH = 1.3;				// used to be 1.1
const double HBONDLENGTH  = 2.5;				// used to be 2.46
const double HBONDANGLE	= 30.0*M_PI/180.0;		// bonding angle has to be less than this value to be considered an H-bond
const double NOBONDLENGTH = 2.0;
const double NHBONDLENGTH = 1.3;		// uhmm... check this?

/* Encoding of the different coordination types
 * The numbering is based on each O having a value of 1, and each H haveing a value of 10 (i.e. add 1 for every O, and 10 for every H...). So a water in a state of OOHH bonding would have a coordination of 22, and a coordination of 13 would be OOOH, 12 = OOH, 11 = OH, 10 = H, etc.
 */
typedef enum {
	UNBOUND=0, O=1, OO=2, OOO=3, OOOO=4, 			// no H
	H=10, OH=11, OOH=12, OOOH=13, OOOOH=14,			// 1 H
	HH=20, OHH=21, OOHH=22, OOOHH=23, OOOOHH=24,		// 2 Hs
	HHH=30, OHHH=31, OOHHH=32, OOOHHH=33, OOOOHHH=34,	// 3 Hs
	HHHH=40, OHHHH=41, OOHHHH=42, OOOHHHH=43, OOOOHHHH=44
} coordination;
// And hopefully that covers all the bonding coordination types :)

typedef std::map<coordination, std::string> coord_map;

// bond types
typedef enum {unbonded, nobond, nhbond, hbond, ohbond, covalent} bondtype;


/********** bond ************/
class Bond {
public:

	Bond ();
	Bond (double length, bondtype btype);
	Bond (const Bond& oldBond);
	~Bond ();

	double	bondlength;	// bond length
	bondtype bond;	 	// bond type (covalent, hydrogen-bond, or unbonded)

	static int num_bonds;

	void SetBondType (const bondtype btype);		// sets the bond type based on the current

};


typedef std::vector<Bond *> Bond_ptr_vec;


#endif
