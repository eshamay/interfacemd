#ifndef CONNECTMATRIX_H_
#define CONNECTMATRIX_H_

#include "utility.h"
#include "bond.h"
#include "atom.h"
#include "h2o.h"

/* The idea of a connectivity matrix is the same as mapping the connections between different members of a system. In this case, we're dealing with atoms in a system, and the connections are defined by the distances between the atoms. We're going to employ a few tricks here in order to make data retrieval a bit more succinct, and also to store more information into one matrix.

Firstly, the matrix is represented by a 2D array of numbers. Each number will represent the distance between 2 atoms. The matrix is symmetric such that each row and each column represent the ordered list of atoms in the system. The upper triangle (all elements above the diagonal) hold the distances between atoms. Only using one half of a symmetric matrix ensures that we don't double our efforts and redo calculations between atoms.

Because we'll be working with water, a hydrogen bond is any bond less than 2.4 angstroms and greater than 1.3 (or so... this is to be defined below). If we want to find the number of hydrogen bonds an atom is involved in, then we look at that atom's column from the top to the diagonal, and from the diagonal to the end of the row, and count the number of bonds that fall within the distance of an H-bond.

As for the diagonal elements - instead of writing a routine that will count the number of h-bonds by looking over the rows and columns, as distances are calculated the diagonal elements for each atom will be updated to reflect the number of H-bonds formed. Thus at any time we can look at the diagonal and know immediately the number of H-bonds an atom is involved in.

This leaves the bottom-diagonal free to store more information. If two atoms are covalently bound, then the bottom diagonal element will mark this with a 1.0.
*/

typedef std::vector<Bond *> Bond_ptr_vec;
typedef std::vector<Bond> Bond_vec;
typedef std::vector< Bond_ptr_vec > Bond_ptr_matrix;
typedef std::vector< Bond_vec > Bond_matrix;
typedef std::vector< std::vector< VecR > > Distance_matrix;

class AdjacencyMatrix {

private:
	Bond_matrix		_matrix;	// the actual data structure for storing connection data
	Atom_ptr_vec	_atoms;		// the atoms in the system
	unsigned int 	_size;
	bool 			_built;		// has the matrix been built to the current size?

public:

	// constructor builds the matrix based on number of atoms to analyze
	AdjacencyMatrix ();
	AdjacencyMatrix (const Atom_ptr_vec& atoms);
	~AdjacencyMatrix ();

	void BuildMatrix ();
	void DeleteMatrix ();
	void UpdateMatrix (const Atom_ptr_vec& atoms, std::string const sys="xyz");
	void _FixSharedAtoms ();
	void ClearBonds ();

	void SetBond (int x, int y, const double blength, const bondtype btype);
	void SetBond (Atom const * const a1, Atom const * const a2, const bondtype btype) {
		Bond * b = GetBond(a1, a2);
		b->bond = btype;
	}

	Bond * GetBond (Atom const * const a1, Atom const * const a2);

	bool CovalentBond (bondtype const b) const;
	int ID (Atom const * const ap) const;

	double Distance (Atom const * const a1, Atom const * const a2) {
		Bond * b = GetBond (a1, a2);
		return (b->bondlength);
	}

	Bond_ptr_vec Bonds (Atom const * const ap);
	Bond_ptr_vec HBonds (Atom const * const ap);
	Atom_ptr_vec BondedAtoms (Atom const * const ap);
	Atom_ptr_vec BondedAtoms (Atom const * const ap, bondtype const b);
	Atom_ptr_vec BondedAtoms (Atom const * const ap, bondtype const b, const string name);

	int NumBonds (Atom const * const ap) {
		return (Bonds(ap).size());
	}
	int NumHBonds (Atom const * const ap) {
		return (HBonds(ap).size());
	}
	int NumHBonds (Water const * const wat);

	coordination WaterCoordination (Water const * const wat);

/*
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

*/
};


#endif
