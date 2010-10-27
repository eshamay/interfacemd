#ifndef GMXSYSTEM_H_
#define GMXSYSTEM_H_

#ifndef CPLUSPLUS
#define CPLUSPLUS 1 
#endif

#include "mdsystem.h"
#include "trrfile.h"
#include "xtcfile.h"
#include "grofile.h"

namespace gromacs {

	//! The template class forms the basis for working with gromacs systems and reading coordinate files. Both TRR and XTC file formats are available. To load either, the appropriate class type has to be used in the template. TRRFile for .trr, and XTCFile for xtc
	template <class T>
		class GMXSystem : public MDSystem {

			public:
				GMXSystem (const std::string gro_filepath, const std::string xdr_filepath);

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

				//! output to stdout some info parsed from the .gro file
				void PrintGROInfo () const { _gro.Print(); }
				//! output (lot's of) info parsed from the xdr/trajectory file
				void PrintXDRInfo () const { _xdr.PrintInfo(); }

			private:
				gromacs::GROFile _gro;
				//! either TRRFile or XTCFile
				T _xdr;		

				void _ParseMolecules ();
		};

	template <class T>
	GMXSystem<T>::GMXSystem (const std::string gro_filepath, const std::string xdrfile_path)
		// ***** need to figure out how to delegate either trr or xtc files here ****
		:
			_gro(gro_filepath), _xdr(xdrfile_path)
	{ 
		this->_ParseMolecules();
	}

	template <class T>
	void GMXSystem<T>::LoadFirst() {
		_xdr.LoadFirst();
		this->_ParseMolecules();
	}

	template <class T>
	void GMXSystem<T>::LoadNext() {
		_xdr.LoadNext();
		this->_ParseMolecules();
	}

	template <>
	void GMXSystem<XTCFile>::_ParseMolecules () {

		VecR_it coord_it = _xdr.begin_coords();

		for (Atom_it it = _gro.begin(); it != _gro.end(); it++) {
			(*it)->Position(*coord_it);
			++coord_it;
		}
	}

	template <>
	void GMXSystem<TRRFile>::_ParseMolecules () {

		VecR_it coord_it = _xdr.begin_coords();
		VecR_it vel_it = _xdr.begin_vels();
		VecR_it force_it = _xdr.begin_forces();

		for (Atom_it it = _gro.begin(); it != _gro.end(); it++) {
			(*it)->Position(*coord_it);
			(*it)->Force(*force_it);

			++coord_it;
			++vel_it;
			++force_it;
		}
	}

} // namespace gromacs

#endif
