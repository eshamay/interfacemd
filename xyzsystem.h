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

	/* some internal methods */
	void _ParseMolecules ();		// take the atoms we have and stick them into molecules
	void _ParseWaters ();
	void _ParseNitrates ();
	void _ParseSulfides ();
	void _ParseWanniers ();
	// Check if an atom pair (O & H) pass the hbond angle criteria
	//bool _HBondAngle (const Atom *H, const Atom *O);
	void _UpdateUnparsedList (Atom_ptr_vec& parsed);
	void _CheckForUnparsedAtoms () const;

public:
	// constructors
	XYZSystem (const std::string& filepath, const VecR& size, const std::string& wannierpath = "");
	~XYZSystem ();

	bondgraph::BondGraph graph;

	// Controller & Calculation methods
	// Update the system to the next timestep
	void LoadNext ();
	void LoadFirst ();
	void Seek (int step);
	int NumSteps () const { return _coords.NumSteps(); }		// number of timesteps in the xyzfile
	int Current () const { return _coords.Current(); }

	void SetReparseLimit (const int limit) { _reparse_limit = limit; }

	const std::vector<VecR>& Wanniers () const { return _wanniers.Coords(); }

	Atom_ptr_vec CovalentBonds (Atom const * const atom) const { return graph.BondedAtoms(atom, bondgraph::covalent); }
	Atom_ptr_vec BondedAtoms (Atom const * const atom) const { return graph.BondedAtoms (atom); }

	VecR SystemDipole ();	// calculate the total system dipole and return it

	typedef std::exception xyzsysex;

	struct unaccountedex : public xyzsysex { 
	  char const* what() const throw() { return "After parsing of the xyz system molecules an atom(s) was left unaccounted for"; }
	};


};	// xyz system class

#endif
