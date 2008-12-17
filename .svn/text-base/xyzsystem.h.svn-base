#ifndef XYZSYSTEM_H_
#define XYZSYSTEM_H_

#include <algorithm>
#include "xyzfile.h"
#include "bondgraph.h"
#include "molecule.h"
#include "h2o.h"
#include "hno3.h"
#include "wannier.h"
#include "utility.h"


class XYZSystem {
	
private:
	XYZFile				_atoms;			// Atomlist parsed from an xyz file
	vector<Molecule *> 	_mols;			// once we have atoms defined, we form them into molecules
	WannierFile 		_wanniers;		// The wannier centers
	BondGraph 			_bondgraph;		// a really useful graph of bonding in the system

	bool _parsed;

	VecR	_dims;			// A vector holding the system size (dimensions)

	//VecR	_centerofmass;	// The total system center of mass

	/* some internal methods */
	void _ParseMolecules ();		// take the atoms we have and stick them into molecules
	void _ParseWanniers ();
	// Check if an atom pair (O & H) pass the hbond angle criteria
	//bool _HBondAngle (const Atom *H, const Atom *O);

public:
	// constructors
	XYZSystem (string filepath, VecR size, string wannierpath);
	~XYZSystem ();
	
	// Controller & Calculation methods
	void LoadNext ();	 					// Update the system to the next timestep 
	void LoadFirst ();
	void Seek (int step);
	int NumSteps () const { return _atoms.NumSteps(); }		// number of timesteps in the xyzfile

	// Input
	//void Dims (VecR size) { _dims = size; }		// Set the system size
	//VecR Dims () const { return _dims; }		// returns the system size. Not the uppercase S to distinguish from the number of atoms

	// Output
	//vector<Atom *>& Atoms () { return _atoms.Atoms(); }	// is this necessary?
	vector<Molecule *>& Molecules () { return _mols; }
	Molecule * Molecules (int mol) { return _mols[mol]; }
	Atom * Atoms (int atom) { 
		Atom * pa;
		pa = _atoms[atom];
		return pa; 
	}
	const vector<VecR>& Wanniers () const { return _wanniers.Coords(); }
	int size ()	const { return _atoms.size(); }

	// returns the distance between two atoms in the system
	double Distance (const Atom * atom1, const Atom * atom2) const { return _bondgraph.Distance (atom1, atom2); }
	
	// calculates (and sets the _centerofmass) the center of mass of the system. For now it only calculates based on waters
	//VecR UpdateCenterOfMass ();
	//VecR CenterOfMass () const { return _centerofmass; }
	//void FindHBonds();	// locates all the H-bonds in a system
	vector<Atom *> CovalentBonds (const Atom * atom) { return _bondgraph.CovalentBonds(atom); }
	std::vector<Atom *> AdjacentAtoms (const Atom * atom) const { return _bondgraph.AdjacentAtoms (atom); }
	std::vector<Atom *> AdjacentAtoms (const Atom * atom, const string name) const { return _bondgraph.AdjacentAtoms (atom, name); }

	VecR SystemDipole ();	// calculate the total system dipole and return it

	int Current ()		const { return _atoms.Current(); }

	// operators
	Atom * operator[] (int index) { 
		Atom * pa;
		pa = _atoms[index];
		return pa;
	}
};

#endif
