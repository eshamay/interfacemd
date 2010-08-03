#pragma once
#ifndef H2O_H_
#define H2O_H_

#include "matrixr.h"
#include "molecule.h"
#include <vector>


// A water class to add a few functions for dealing with water molecules specifically.
class Water: public Molecule {

protected:
	//double _eulerAngles[2][3];
	VecR _oh1, _oh2;				// Both of the OH vectors

	AtomPtr _o, _h1, _h2;			// pointers to the atoms for easy access

public:
	Water ();	// a default constructor
	~Water ();	// a destructor
	Water (const Molecule& molecule);	// copy constructor for casting from a molecule
	Water (const MolPtr& mol);

	typedef Water* WaterPtr;
	static int numWaters;			// total number of waters in the system

	// Functions for analysis
	void SetAtoms ();
	void SetBondLengths ();

	VecR Bisector ();		// calculates the bisector (unit vector) of the water
	VecR Normal () const { return _y; }

	// flip the molecule about a plane (i.e. take its mirror image about a plane) through the oxygen about a given axis
	void Flip (const coord axis);

	void SetMoritaAxes (const int bond = 1);		// Determines the molecular-frame axes (a la Morita&Hynes2000) with one bond on the Z-axis, the other in the positive X direction.
	void SetOrderAxes ();		// sets the molecular axes such that the z-axis is along the C2V axis point from the H's to the O, and the x-axis is in the plane of the molecule

	MatR const & DCMToLab ();
	MatR const & DCMToLabMorita (const int bond = 1);	// get the direction cosine matrix for rotations to the lab frame from the morita-hynes one
	MatR const & DCMToLabOrder ();						// direction cosine matrix using the bisector as the molecular z-axis

	MatR EulerMatrix;					// The euler rotation matrix
	double EulerAngles[3];				// euler angles as defined in "The Raman Effect" Appendix A5 (theta, phi, chi)

	void CalcEulerAngles ();

	#ifdef WATER_POLARIZ
	void CalcAlpha ();		// calculate the molecular polarizability tensor (as per morita+hynes 2002 method)
	MatR const & Alpha () const { return _alpha; }
	#endif

	AtomPtr O () { return _o; }
	AtomPtr H1 () { return _h1; }
	AtomPtr H2 () { return _h2; }

	VecR const * OH1 () const { return &_oh1; }
	VecR const * OH2 () const { return &_oh2; }
	double Angle () const { return (_oh1 < _oh2); }	// returns the cos of the H-O-H angle

	VecR MolecularAxis (); 
};

typedef Water::WaterPtr WaterPtr;
typedef std::vector<WaterPtr> Water_ptr_vec;
typedef std::vector<Water> Water_vec;


#endif
