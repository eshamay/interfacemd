#ifndef AMBERSYSTEM_H_
#define AMBERSYSTEM_H_

#include "mdsystem.h"
#include "crdfile.h"
#include "forcefile.h"
#include "topfile.h"

class AmberSystem : public MDSystem {

private:
	Atom_ptr_vec	_atoms;
	TOPFile		_topfile;
	CRDFile		_coords;
	ForceFile	_forces;

	void _ParseAtomInformation ();
	void _ParseAtomVectors ();
	void _ParseMolecules ();

public:
	// constructors
	AmberSystem () { }
	AmberSystem (const std::string prmtop, const std::string mdcrd, const std::string mdvel);
	~AmberSystem ();

	// Controller & Calculation methods
	void LoadNext ();	 					// Update the system to the next timestep
	void LoadFirst ();

	bool eof () { return _coords.eof(); }

	// Output
	VecR	Dims () 		const 	{ return _coords.Dims(); }		// returns the system size.

	int 	Current ()		const 	{ return _coords.Current(); }

	void PrintCRDFile () const;						// to output a frame of the system in .crd format

	VecR& Forces (int index)	{ return _forces[index]; }
	VecR& Coords (int index)	{ return _coords[index]; }

};

#endif
