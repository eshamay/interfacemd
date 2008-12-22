#ifndef H2O_H_
#define H2O_H_

#include "molecule.h"

#ifdef H2O_DIPOLE_PARM
#include "dipoleparm.h"
#endif

#include "matrixr.h"
#include "utility.h"


#define R_eq	0.9575		// equilibrium OH bond length
#define Theta_eq	104.51	// equilibrium H-O-H angle
#define qO_eq	-0.6750		// equilibrium partial charges of the oxygen
#define qH_eq	 0.3285		// 		and the hydrogens


// A water class to add a few functions for dealing with water molecules specifically.
class Water: public Molecule {

protected:
	//double _eulerAngles[2][3];
	VecR _oh1, _oh2;				// Both of the OH vectors

	Atom *_o, *_h1, *_h2;			// pointers to the atoms for easy access
	#ifdef WATER_POLARIZ
	MatR _alpha;					// polarizability of the molecule
	#endif

#ifdef H2O_DIPOLE_PARM
	static WaterDipoleParms _dipparms;		// The water dipole parameter file
#endif

public:
	Water ();	// a default constructor
	~Water ();	// a destructor
	Water (const Molecule& molecule);	// copy constructor for casting from a molecule

	static int numWaters;			// total number of waters in the system
	// Functions for analysis
	void SetAtoms ();
	VecR Bisector ();		// calculates the bisector (unit vector) of the water
	VecR Normal () const { return _y; }

	void CalcDipole ();
	VecR const & Dipole () const { return _dipole; }			// calculates the dipole (from a parameterized source)
	
	void SetMoritaAxes (const int bond = 1);		// Determines the molecular-frame axes (a la Morita&Hynes2000) with one bond on the Z-axis, the other in the positive X direction.
	void SetOrderAxes ();

	MatR const & DCMToLab ();							// get the direction cosine matrix for rotations to the lab frame
	MatR const & DCMToLabMorita (const int bond = 1);	// get the direction cosine matrix for rotations to the lab frame from the morita-hynes one
	MatR const & DCMToLabOrder ();						// direction cosine matrix using the bisector as the molecular z-axis
	MatR DCM;											// the direction cosine matrix for rotating the molecule to the lab frame

	MatR EulerMatrix;					// The euler rotation matrix
	double EulerAngles[3];				// euler angles as defined in "The Raman Effect" Appendix A5 (theta, phi, chi)

	void CalcEulerAngles ();

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
