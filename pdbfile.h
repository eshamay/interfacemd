#ifndef PDBFILE_H_
#define PDBFILE_H_

#include <string>
#include <vector>
#include <iostream>
#include "mdsystem.h"

using namespace std;

class PDBFile : public MDSystem {

	FILE *_file;				// the PDB file listing all the atom coordinates
	string _path;

	int _laststep,
		_currentstep;

	int	_numAtoms;				// total number of atoms in the system
	int _numMols;

	int _loaded;				// To tell wether or not a file has been loaded

	int _FindLastStep();		// run through the file and find the last frame available
	Atom *_ParseAtom (const char *line);
	void _ParseMolecules();

public:

	PDBFile (string path);
	PDBFile (std::vector<Molecule *>& mols);
	PDBFile ();
	~PDBFile ();

	// Various control functions
	void LoadNext ();
	void LoadFirst ();
	void LoadLast ();
	void Seek (int step);

	// output functions
	int Last () { return _laststep; }
	int Current () { return _currentstep; }

	static void WritePDB (std::vector<Molecule *>& system);		// given a vector of molecules, this will print out a PDB file

};

#endif
