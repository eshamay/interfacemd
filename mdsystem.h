#pragma once
#ifndef MDSYSTEM_H_
#define MDSYSTEM_H_

#include "vecr.h"
#include "atom.h"
#include "molecule.h"
#include <math.h>
#include <string>
#include <vector>

class MDSystem {

	protected:
		Atom_ptr_vec	_atoms;		// the atoms in the system
		Mol_ptr_vec		_mols;		// the molecules in the system

		static VecR		_dimensions;		// system dimensions - size
		
		//! Parses out molecules from the set of atoms in an MD data set. This is typically done via topology files, or some other defined routine that determines connectivity between atoms to form molecules.
		virtual void _ParseMolecules () = 0;
		bool _parse_molecules;	// this gets set if the molecules are to be parsed to determine the specific types.

	public:

		virtual ~MDSystem();

		//! Rewinds the necessary input and data files and then loads the first frame of an MD simulation
		virtual void LoadFirst () = 0;
		//! Loads the next frame of an MD simulation data set
		virtual void LoadNext () = 0;


		//! The set of all molecules in a system
		Mol_ptr_vec& Molecules () { return _mols; }
		//! An iterator to the beginning of the set of molecules
		Mol_it begin_mols () const { return _mols.begin(); }
		//! An iterator to the end of the set of molecules
		Mol_it end_mols () const { return _mols.end(); }
		//! An indexing method for retrieving specific molecules in a system
		MolPtr Molecules (int index) { return _mols[index]; }
		//! Returns the total number of molecules in a system
		int NumMols () const { return _mols.size(); }

		Atom_ptr_vec& Atoms () { return _atoms; }
		Atom_it begin () { return _atoms.begin(); }
		Atom_it end () { return _atoms.end(); }
		AtomPtr Atoms (const int index) { return _atoms[index]; }
		AtomPtr operator[] (int index) { return _atoms[index]; }
		int NumAtoms ()	const { return (int)_atoms.size(); }

		int size () const { return (int)_atoms.size(); }
		static VecR Dimensions () { return MDSystem::_dimensions; }
		static void Dimensions (const VecR& dimensions) { MDSystem::_dimensions = dimensions; }

		/* Beyond simple system stats, various computations are done routinely in a molecular dynamics system: */

		// Calculate the distance between two points within a system that has periodic boundaries
		static VecR Distance (const VecR v1, const VecR v2) ;

		// Calculate the distance between two atoms given the periodic boundaries of the system
		static VecR Distance (const AtomPtr atom1, const AtomPtr atom2) ;

		// Calculates the minimum distance between two molecules - i.e. the shortest inter-molecular atom-pair distance
		static double Distance (const MolPtr mol1, const MolPtr mol2) ;

		//! Calculates a molecular dipole moment using the "classical E&M" method - consider each atom in the molecule as a point-charge, and that the molecule has no net charge (calculation is independent of origin location). This returns a vector that is the sum of the position*charge (r*q) of each atom.
		static VecR CalcClassicDipole (MolPtr mol);

		//! Calculates the dipole of a molecule using the classic E&M method, but also account for the wannier localization centers (WLC). This assumes that the wannier localization centers have also been calculated and parsed into the molecule.
		static VecR CalcWannierDipole (MolPtr mol);
};

#include "moleculefactory.h"

#endif
