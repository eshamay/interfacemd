#ifndef XYZSYSTEM_H_
#define XYZSYSTEM_H_

#include <algorithm>
#include "mdsystem.h"
#include "xyzfile.h"
#include "wannier.h"
#include "graph.h"

#define DEBUG 1

class XYZSystem : public MDSystem {

private:
	XYZFile				_coords;		// Atomlist parsed from an xyz file
	WannierFile 		_wanniers;		// The wannier centers
	bool _parsed;
	BondGraph		_graph;			// a really useful graph for playing with atoms and bonds

	/* For debugging (and other useful things?) this will keep a list of all the atoms that have been processed into molecules. Any atoms left over at the end of the parsing routine are not included and ... can potentially cause problems */
	Atom_ptr_vec _unparsed;

	/* some internal methods */
	void _ParseMolecules ();		// take the atoms we have and stick them into molecules
	void _ParseWaters ();
	void _ParseNitrates ();
	void _ParseWanniers ();
	// Check if an atom pair (O & H) pass the hbond angle criteria
	//bool _HBondAngle (const Atom *H, const Atom *O);

public:
	// constructors
	XYZSystem () { _graph.SysType("xyz"); }
	XYZSystem (string filepath, VecR size, string wannierpath);
	~XYZSystem ();

	// Controller & Calculation methods
	// Update the system to the next timestep
	void LoadNext ();
	void LoadFirst ();
	void Seek (int step);
	int NumSteps () const { return _coords.NumSteps(); }		// number of timesteps in the xyzfile
	int Current ()		const { return _coords.Current(); }

	const std::vector<VecR>& Wanniers () const { return _wanniers.Coords(); }

	Atom_ptr_vec CovalentBonds (Atom const * const atom) const { return _graph.BondedAtoms(atom, covalent); }
	Atom_ptr_vec BondedAtoms (Atom const * const atom) const { return _graph.BondedAtoms (atom); }

	VecR SystemDipole ();	// calculate the total system dipole and return it
};

#endif
