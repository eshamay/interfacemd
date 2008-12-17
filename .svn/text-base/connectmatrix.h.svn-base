#ifndef CONNECTMATRIX_H_
#define CONNECTMATRIX_H_

#include "h2o.h"
#include "utility.h"
#include <algorithm>


/* The idea of a connectivity matrix is the same as mapping the connections between different members of a system. In this case, we're dealing with atoms in a system, and the connections are defined by the distances between the atoms. We're going to employ a few tricks here in order to make data retrieval a bit more succinct, and also to store more information into one matrix.

Firstly, the matrix is represented by a 2D array of numbers. Each number will represent the distance between 2 atoms. The matrix is symmetric such that each row and each column represent the ordered list of atoms in the system. The upper triangle (all elements above the diagonal) hold the distances between atoms. Only using one half of a symmetric matrix ensures that we don't double our efforts and redo calculations between atoms.

Because we'll be working with water, a hydrogen bond is any bond less than 2.4 angstroms and greater than 1.3 (or so... this is to be defined below). If we want to find the number of hydrogen bonds an atom is involved in, then we look at that atom's column from the top to the diagonal, and from the diagonal to the end of the row, and count the number of bonds that fall within the distance of an H-bond.

As for the diagonal elements - instead of writing a routine that will count the number of h-bonds by looking over the rows and columns, as distances are calculated the diagonal elements for each atom will be updated to reflect the number of H-bonds formed. Thus at any time we can look at the diagonal and know immediately the number of H-bonds an atom is involved in.

This leaves the bottom-diagonal free to store more information. If two atoms are covalently bound, then the bottom diagonal element will mark this with a 1.0. 
*/

const double OH_BOND_LENGTH = 1.3;
const double H_BOND_LENGTH  = 2.46;
const double H_BOND_ANGLE	= 0.866025;		// cos(theta) has to be less than this value to be considered an H-bond
const double NO_BOND_LENGTH = 2.0;
const double NH_BOND_LENGTH = 1.3;		// uhmm... check this?

//For naming the bond-types, let's initialize an array with readable coordination names;
const string COORD_NAMES[14] = {"Unbound", "O", "OH", "OHH", "OO", "OOH", "OOHH", "H", "HH", "OOOH", "OOOHH", "OOHHH", "OVER", "OTHER"};

// Encoding of the different coordination types
enum coordination {UNBOUND=0, O, OH, OHH, OO, OOH, OOHH, H, HH, OOOH, OOOHH, OOHHH, OVER, OTHER};

// A non-member function for getting a nice string representation of a coordination
const string CoordName (coordination coord);

class ConnectivityMatrix {

private:
	std::vector< std::vector<double> >	_matrix;		// the actual data structure for storing connection data
	std::vector<Atom *>				_atoms;			// the atoms in the system

	// forms a covalent bond (bottom-triangle of the matrix) between two atoms
	void _FormCovalentBond (const Atom * atom1, const Atom * atom2);

	// let's us update the H-bond information (diagonal elements) from abroad
	void _FormHBond (Atom * atom1, Atom * atom2);

	// runs through bonding criteria of different atom pairs
	void _FindBonds (Atom * atom1, Atom * atom2);

public:
	// constructor
	ConnectivityMatrix (std::vector<Atom *>& atoms);

	void UpdateMatrix ();

	// returns the number of H-bonds that an atom is involved in (i.e. the diagonal element)
	int HBonds (const Atom * atom) const { 
		return (int)(_matrix[atom->ID()][atom->ID()]); 
	}	
	
	coordination FindWaterCoordination (const Water& water) const;	// returns the coordination number (defined above) of a water

	// return a list of all the atoms that are covalently bound to a given atom
	std::vector<Atom *> CovalentBonds (const Atom * atom) const;
	
	// returns the distance between two atoms
	double Distance (const Atom * atom1, const Atom * atom2) const;		
	int    size() const { return _atoms.size(); }

	// returns the closest atoms of a given name to a given atom
	// input is the target atom's id, atomname is the name of the other atoms in the system we want returned,
	// and number is the number of nearest atoms 
	// i.e. ClosestAtoms (5, O, 3) - returns the three closest O's to the atom with ID 5
	std::vector<Atom *> ClosestAtoms (const int input, const string atomname, const int number) const;

};
	

#endif
