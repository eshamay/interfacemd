#ifndef MDSYSTEM_H_
#define MDSYSTEM_H_

#include <vector>
#include <math.h>
#include <string>
#include "vecr.h"
#include "atom.h"
#include "molecule.h"
#include "utility.h"

typedef std::vector<Atom *> Atom_ptr_vec;
typedef std::vector<Molecule *> Mol_ptr_vec;

class MDSystem {

protected:
	Atom_ptr_vec	_atoms;		// the atoms in the system
	Mol_ptr_vec		_mols;		// the molecules in the system

	static VecR		_dimensions;		// system dimensions - size

public:

	virtual ~MDSystem();

	virtual void LoadNext () = 0;
	virtual void LoadFirst () = 0;

	virtual void _ParseMolecules () = 0;

	Mol_ptr_vec& Molecules () { return _mols; }
	Molecule * Molecules (int index) { return _mols[index]; }
	int NumMols () const { return _mols.size(); }

	Atom_ptr_vec& Atoms () { return _atoms; }
	Atom * Atoms (const int index) { return _atoms[index]; }
	Atom * operator[] (int index) { return _atoms[index]; }
	int NumAtoms ()	const { return (int)_atoms.size(); }

	int size () const { return (int)_atoms.size(); }
	static VecR Dimensions () { return _dimensions; }
	static void Dimensions (const VecR& dimensions) { MDSystem::_dimensions = dimensions; }

	/* Beyond simple system stats, various computations are done routinely in a molecular dynamics system: */

	// Calculate the distance between two points within a system that has periodic boundaries
	static VecR Distance (const VecR& v1, const VecR& v2);

	// Calculate the distance between two atoms given the periodic boundaries of the system
	static VecR Distance (const Atom * atom1, const Atom * atom2);
};

#include "h2o.h"
#include "h3o.h"
#include "hno3.h"
#include "carbonchain.h"
#include "decane.h"
#include "oh.h"

#endif
