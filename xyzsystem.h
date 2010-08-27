#pragma once
#ifndef XYZSYSTEM_H_
#define XYZSYSTEM_H_

#include "mdsystem.h"
#include "xyzfile.h"
#include "wannier.h"
#include "graph.h"
#include <ext/algorithm>
#include <exception>

class XYZSystem : public MDSystem {

	private:
		XYZFile				_coords;		// Atomlist parsed from an xyz file
		WannierFile 		_wanniers;		// The wannier centers

		// This is set to control how often the system molecules will be reparsed.
		// Note that the wanniers centers need to be loaded with every timestep,
		// but the molecules don't necessarily need to
		int _reparse_limit;					
		int _reparse_step;


		/* For debugging (and other useful things?) this will keep a list of all the atoms that have been processed into molecules. Any atoms left over at the end of the parsing routine are not included and ... can potentially cause problems */
		Atom_ptr_vec _unparsed;

		void _ParseMolecules ();		// take the atoms we have and stick them into molecules - general umbrella routine

		// Parses a simple molecule that is composed of a central atom, and all other atoms are connected to it - i.e. h2o, no3, ccl4, etc
		template <typename T>
			void _ParseSimpleMolecule (const Atom::Element_t central_elmt, const Atom::Element_t outer_elmt, const unsigned numOuter);

		/* mor specialized parsing routines */

		void _ParseNitricAcids ();
		void _ParseProtons ();
		void _ParseWanniers ();


		void _UpdateUnparsedList (Atom_ptr_vec& parsed);	// fixes the list of unparsed atoms
		void _CheckForUnparsedAtoms () const;

	public:
		// constructors
		XYZSystem (const std::string& filepath, const VecR& size, const std::string& wannierpath = "");
		~XYZSystem ();

		bondgraph::BondGraph graph;

		// Update the system to the next timestep
		void LoadNext ();
		void LoadFirst ();
		void Seek (int step);
		int NumSteps () const { return _coords.NumSteps(); }		// number of timesteps in the xyzfile
		int Current () const { return _coords.Current(); }

		void SetReparseLimit (const int limit) { _reparse_limit = limit; }

		const VecR_vec& Wanniers () const { return _wanniers.Coords(); }

		Atom_ptr_vec CovalentBonds (const AtomPtr atom) const { return graph.BondedAtoms(atom, bondgraph::covalent); }
		Atom_ptr_vec BondedAtoms (const AtomPtr atom) const { return graph.BondedAtoms (atom); }

		VecR SystemDipole ();	// calculate the total system dipole and return it

		typedef std::exception xyzsysex;

		struct unaccountedex : public xyzsysex { 
			char const* what() const throw() { return "After parsing of the xyz system molecules an atom(s) was left unaccounted for"; }
		};

		template <typename U>
			class AtomPtr_In_List : public std::binary_function<AtomPtr,U,bool> {
				public:
					bool operator () (const AtomPtr ap, const U u) const {
						return find(u.begin(), u.end(), ap) != u.end();
					}
			};



};	// xyz system class



template <typename T>
void XYZSystem::_ParseSimpleMolecule (const Atom::Element_t central_elmt, const Atom::Element_t outer_elmt, const unsigned numOuter) {

	for (Atom_it it = _atoms.begin(); it != _atoms.end(); it++) {

		if ((*it)->Element() != central_elmt) continue;

		// for every S in the system, see if 2 oxygens are connected
		Atom_ptr_vec outers = graph.BondedAtoms (*it, bondgraph::covalent, outer_elmt);
		if (outers.size() != numOuter) continue;

		int molIndex = _mols.size();

		outers.push_back(*it);
		T * newmol = new T();

		_mols.push_back (newmol);
		newmol->MolID (molIndex);

		for (Atom_it jt = outers.begin(); jt != outers.end(); jt++) {
			newmol->AddAtom (*jt);
		}

		_UpdateUnparsedList(outers);
	}

	return;
}	// parse simple molecule	


#endif
