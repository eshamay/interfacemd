#ifndef MDSYSTEM_H_
#define MDSYSTEM_H_

#include <vector>
#include <math.h>
#include <string>
#include "molecule.h"
#include "atom.h"
#include "h2o.h"
#include "h3o.h"
#include "hno3.h"
#include "oh.h"
#include "utility.h"

typedef std::vector<Atom *> Atom_ptr_vec;
typedef std::vector<Molecule *> Mol_ptr_vec;

class MDSystem {

protected:
	Atom_ptr_vec	_atoms;		// the atoms in the system
	Mol_ptr_vec		_mols;		// the molecules in the system

	VecR			_dims;		// system dimensions - size
	int				_numTimeSteps;		// number of timesteps in the MD data files

public:

	virtual void LoadNext () = 0;
	virtual void LoadFirst () = 0;

	virtual void _ParseMolecules () = 0;

	Mol_ptr_vec& Molecules () { return _mols; }
	Molecule * Molecules (int index) { return _mols[index]; }
	const int NumMols () const { return _mols.size(); }

	Atom_ptr_vec& Atoms () { return _atoms; }
	Atom * Atoms (const int index) { return _atoms[index]; }
	Atom * operator[] (int index) { return _atoms[index]; }
	const int NumAtoms ()	const { return _atoms.size(); }

	const int size () const { return (int)_atoms.size(); }
	virtual const int NumSteps () const { return _numTimeSteps; }

};

#endif
