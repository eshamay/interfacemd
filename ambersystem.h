#ifndef AMBERSYSTEM_H_
#define AMBERSYSTEM_H_

#include <vector>
#include <math.h>
#include <string>
#include "crdfile.h"
#include "forcefile.h"
#include "topfile.h"
#include "atom.h"
#include "molecule.h"
#include "h2o.h"
#include "hno3.h"
#include "bondgraph.h"
	
class AmberSystem {

private:	
	TOPFile		_topfile;
	CRDFile		_coords;
	ForceFile	_forces;

	std::vector<Atom *> _atoms;			// A listing of all the atoms in the system with information parsed from the topfile and crdfile
	std::vector<Molecule *> _mols;		// The molecules within a system - defined by the residues in the topology files

	BondGraph	_graph;					// our adjacency graph for calculating connections in the system

	void _ParseAtomInformation ();
	void _ParseAtomVectors ();
	void _ParseMolecules ();

public:
	// constructors
	AmberSystem (string prmtop, string mdcrd, string mdvel);
	~AmberSystem ();
	
	// Controller & Calculation methods
	void LoadNext ();	 					// Update the system to the next timestep 
	void LoadFirst ();

	bool eof () { return _coords.eof(); }

	// Output
	int	 	NumAtoms ()		const 	{ return _atoms.size(); }		// return the number of atoms.
	int		size ()			const 	{ return _atoms.size(); }
	int		NumMols ()		const 	{ return _mols.size(); }		// returns the number of residues
	VecR	Dims () 		const 	{ return _coords.Dims(); }		// returns the system size.

	int 	Current ()		const 	{ return _coords.Current(); }

	void PrintCRDFile (string filepath);						// to output a frame of the system in .crd format

	Molecule * Molecules (int index) { return _mols[index]; }
	std::vector<Molecule *>& Molecules () { return _mols; }
	Atom * Atoms (int index) 	{ return _atoms[index]; }
	std::vector<Atom *>& Atoms () { return _atoms; }
	VecR& Forces (int index)	{ return _forces[index]; }
	VecR& Coords (int index)	{ return _coords[index]; }

	// operators
	Atom * operator[] (int index) { return _atoms[index]; }

	void UpdateGraph (vector<Atom *>& int_atoms) { _graph.UpdateGraph(int_atoms); }						// update the connectivity matrix/bond graph
	//double Distance (Atom * atom1, Atom * atom2) const { return _graph.Distance (atom1, atom2); }
	//vector<Atom *> AdjacentAtoms (Atom * atom) const { return _graph.AdjacentAtoms (atom); }
	//coordination WaterCoordination (Water * wat) const { return _graph.FindWaterCoordination (wat); }
	
	// get the number of H-bonds the atom is involved in
	//int HBonds (const int index) const { return _matrix->HBonds(_atoms[index]); }	

	//coordination WaterCoord (const Water& water) const { return _matrix->FindWaterCoordination(water); }

};

#endif
