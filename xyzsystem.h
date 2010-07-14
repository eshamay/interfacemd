#ifndef XYZSYSTEM_H_
#define XYZSYSTEM_H_

#include <ext/algorithm>
#include "mdsystem.h"
#include "xyzfile.h"
#include "wannier.h"
#include "graph.h"

class XYZSystem : public MDSystem {

private:
	XYZFile				_coords;		// Atomlist parsed from an xyz file
	WannierFile 		_wanniers;		// The wannier centers
	bool _parsed;


	/* For debugging (and other useful things?) this will keep a list of all the atoms that have been processed into molecules. Any atoms left over at the end of the parsing routine are not included and ... can potentially cause problems */
	Atom_ptr_vec _unparsed;

	/* some internal methods */
	void _ParseMolecules ();		// take the atoms we have and stick them into molecules
	void _ParseWaters ();
	void _ParseNitrates ();
	void _ParseSulfides ();
	void _ParseWanniers ();
	// Check if an atom pair (O & H) pass the hbond angle criteria
	//bool _HBondAngle (const Atom *H, const Atom *O);
	void _UpdateUnparsedList (const Atom_ptr_vec& parsed);
	void _CheckForUnparsedAtoms () const;

public:
	// constructors
	XYZSystem (std::string filepath, VecR size, std::string wannierpath = "");
	~XYZSystem ();

	BondGraph graph;

	// Controller & Calculation methods
	// Update the system to the next timestep
	void LoadNext ();
	void LoadFirst ();
	void Seek (int step);
	int NumSteps () const { return _coords.NumSteps(); }		// number of timesteps in the xyzfile
	int Current () const { return _coords.Current(); }

	const std::vector<VecR>& Wanniers () const { return _wanniers.Coords(); }

	Atom_ptr_vec CovalentBonds (Atom const * const atom) const { return graph.BondedAtoms(atom, covalent); }
	Atom_ptr_vec BondedAtoms (Atom const * const atom) const { return graph.BondedAtoms (atom); }

	VecR SystemDipole ();	// calculate the total system dipole and return it

};

#endif
