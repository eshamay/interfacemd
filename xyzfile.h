#ifndef XYZFILE_H_
#define XYZFILE_H_

#include <vector>
#include <string>
#include <stdlib.h>
#include "atom.h"
#include "vecr.h"
#include "utility.h"


class XYZFile {

	std::vector<Atom *> _atoms;		// The listing of the atoms in the file

	FILE *_file;				// the XYZ file listing all the atom coordinates
	std::string _path;

	int _currentstep, _firstStep, _lastStep, _numSteps,
		_numatoms;				// total number of centers to process from the file for the frame

	bool _initialized;				// To tell wether or not a file has been loaded

	void _FindSteps ();			// assuming that each timestep is headed by an "i = ..." line, then we can load info on the first, last, and total timesteps

	/* Internal class methods */
	//Atom _ParseAtom (Atom * pAtom);

public:

	XYZFile (std::string path);
	XYZFile ();
	~XYZFile ();

	// Various control functions
	// see LoadFirst for the init arg
	void LoadNext ();

	// If the atomlist has been initialized already, then send true, otherwise to form a new list, send false.
	void LoadFirst ();
	void Seek (int step);

	// output functions
	std::vector<Atom *>& Atoms () { return _atoms; }
	int Current () const { return _currentstep; }
	int NumSteps () const { return _numSteps; }
	int size () const { return _atoms.size(); }

	Atom * operator[] (int index) { return _atoms[index]; }

};

#endif
