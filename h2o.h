#ifndef H2O_H_
#define H2O_H_

#include "molecule.h"
#include <map>

#ifdef H2O_DIPOLE_PARM
#include "dipoleparm.h"
#endif

#include "matrixr.h"
#include "utility.h"

#define R_eq	0.9575		// equilibrium OH bond length
#define Theta_eq	104.51	// equilibrium H-O-H angle
#define qO_eq	-0.6750		// equilibrium partial charges of the oxygen
#define qH_eq	 0.3285		// 		and the hydrogens

/* Encoding of the different coordination types
 * The numbering is based on each O having a value of 1, and each H haveing a value of 10 (i.e. add 1 for every O, and 10 for every H...). So a water in a state of OOHH bonding would have a coordination of 22, and a coordination of 13 would be OOOH, 12 = OOH, 11 = OH, 10 = H, etc.
 */
typedef enum {
	UNBOUND=0, 
	O=1, 
	OO=2, 
	OOO=3, 
	H=10, 
	OH=11, 
	OOH=12,
	OOOH=13,
	HH=20,
	OHH=21,
	OOHH=22,
	OOOHH=23,
	HHH=30,
	OHHH=31,
	OOHHH=32,
	OOOHHH=33
} coordination;
// And hopefully that covers all the bonding coordination types :)

// A water class to add a few functions for dealing with water molecules specifically.
class Water: public Molecule {

protected:
	//double _eulerAngles[2][3];
	VecR _oh1, _oh2;				// Both of the OH vectors

	Atom *_o, *_h1, *_h2;			// pointers to the atoms for easy access
	#ifdef WATER_POLARIZ
	MatR _alpha;					// polarizability of the molecule
	#endif

	coordination _coord;				// the bonding coordation of the water

	#ifdef H2O_DIPOLE_PARM
	static WaterDipoleParms _dipparms;		// The water dipole parameter file
	#endif

public:
	Water ();	// a default constructor
	~Water ();	// a destructor
	Water (const Molecule& molecule);	// copy constructor for casting from a molecule

	static int numWaters;			// total number of waters in the system
	static map<coordination, string> CoordinationNames;		// a map for printing out names of the different types of coordinations

	// Functions for analysis
	void SetAtoms ();
	VecR Bisector ();		// calculates the bisector (unit vector) of the water
	VecR Normal () const { return _y; }

	void CalcDipole ();
	VecR const & Dipole () const { return _dipole; }			// calculates the dipole (from a parameterized source)
	
	void SetMoritaAxes (const int bond = 1);		// Determines the molecular-frame axes (a la Morita&Hynes2000) with one bond on the Z-axis, the other in the positive X direction.
	void SetOrderAxes ();		// sets the molecular axes such that the z-axis is along the C2V axis point from the H's to the O, and the x-axis is in the plane of the molecule

	MatR const & DCMToLab (const coord axis = z);							// get the direction cosine matrix for rotations to the lab frame
	MatR const & DCMToLabMorita (const coord axis = z);	// get the direction cosine matrix for rotations to the lab frame from the morita-hynes one
	MatR const & DCMToLabOrder ();						// direction cosine matrix using the bisector as the molecular z-axis
	MatR DCM;											// the direction cosine matrix for rotating the molecule to the lab frame

	MatR EulerMatrix;					// The euler rotation matrix
	double EulerAngles[3];				// euler angles as defined in "The Raman Effect" Appendix A5 (theta, phi, chi)

	void CalcEulerAngles (const coord axis = z);

	// dealing with the bonding coordination of a water molecule - how is it hydrogen-bonded?
	void Coordination (const coordination coord) { _coord = coord; }
	coordination Coordination () const { return _coord; }

	#ifdef WATER_POLARIZ
	void CalcAlpha ();		// calculate the molecular polarizability tensor (as per morita+hynes 2002 method)
	MatR const & Alpha () const { return _alpha; }
	#endif
	
	VecR const * OH1 () const { return &_oh1; }
	VecR const * OH2 () const { return &_oh2; }
	double Angle () const { return acos(_oh1 < _oh2) * 180.0/M_PI; }
	
};


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
	VecR const * OH1 () const { return &_oh1; }
	VecR const * OH2 () const { return &_oh2; }
	VecR const * OH3 () const { return &_oh3; }
};


#endif
