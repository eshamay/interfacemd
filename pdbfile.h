#ifndef PDBFILE_H_
#define PDBFILE_H_

#include <string>
#include <vector>
#include <iostream>
#include "atom.h"
#include "molecule.h"
#include "utility.h"

using namespace std;

class PDBFile {

	vector<Molecule *> _molecules;		// The listing of the molecules in the file
	vector<Atom *> _atoms;

	FILE *_file;				// the PDB file listing all the atom coordinates
	string _path;

	int _initstep,				// the first frame listed in the file
		_laststep,
		_currentstep;

	int	_numAtoms;				// total number of atoms in the system
	int _numMols;

	int _loaded;				// To tell wether or not a file has been loaded

	int _FindLastStep();		// run through the file and find the last frame available
	Atom *_ParseAtom (const char *line);

public:

	PDBFile (string path);
	PDBFile (vector<Molecule *>& mols);
	PDBFile ();
	~PDBFile ();

	// Various control functions
	void LoadNext ();
	void LoadFirst ();
	void LoadLast ();
	void Seek (int step);

	// output functions
	int First () { return _initstep; }
	int Last () { return _laststep; }
	int Current () { return _currentstep; }
	int size () { return _numAtoms; }
	int numAtoms () { return _numAtoms; }
	int numMols () { return _numMols; }

	static void WritePDB (vector<Molecule *>& system);		// given a vector of molecules, this will print out a PDB file

	vector<Molecule *>& Molecules () { return _molecules; }
	Atom * Atoms (int index) { return _atoms[index]; }
	Molecule *operator[] (int index) { return _molecules[index]; }

};

#endif
