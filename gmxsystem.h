#ifndef GMXSYSTEM_H_
#define GMXSYSTEM_H_

#ifndef CPLUSPLUS
#define CPLUSPLUS 1 
#endif

#include "mdsystem.h"
#include "trrfile.h"
#include "grofile.h"


class GMXSystem : public MDSystem {

	public:
		GMXSystem (const std::string gro_filepath, const std::string trr_filepath);

		void LoadNext ();
		void LoadFirst ();

		Atom_ptr_vec& Atoms () { return _gro.Atoms(); }
		AtomPtr Atoms (const int index) { return _gro.Atoms(index); }
		Atom_it begin () { return _gro.begin(); }
		Atom_it end () { return _gro.end(); }

		Mol_ptr_vec& Molecules () { return _gro.Molecules(); }
		Mol_it begin_mols () const { return _gro.begin_mols(); }
		Mol_it end_mols () const { return _gro.end_mols(); }
		MolPtr Molecules (int index) { return _gro.Molecules(index); }
		int NumMols () const { return _gro.NumMols(); }

		AtomPtr operator[] (int index) { return _gro[index]; }
		int NumAtoms ()	const { return (int)_gro.NumAtoms(); }
		int size () const { return (int)_gro.size(); }

	private:
		gromacs::GROFile _gro;
		gromacs::TRRFile _trr;

		void _ParseMolecules ();
};

#endif
