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
	virtual ~Water ();	// a destructor
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

	AtomPtr O () { return _o; }
	AtomPtr H1 () { return _h1; }
	AtomPtr H2 () { return _h2; }

	virtual VecR ReferencePoint () const { return _o->Position(); }

	VecR const * OH1 () const { return &_oh1; }
	VecR const * OH2 () const { return &_oh2; }
	double Angle () const { return (_oh1 < _oh2); }	// returns the cos of the H-O-H angle

	virtual void SetOrderAxes ();		// sets the molecular axes such that the z-axis is along the C2V axis point from the H's to the O, and the x-axis is in the plane of the molecule
	virtual VecR MolecularAxis (); 
};

typedef Water::WaterPtr WaterPtr;
typedef std::vector<WaterPtr> Water_ptr_vec;
typedef Water_ptr_vec::iterator Wat_it;
typedef std::vector<Water> Water_vec;


#endif
