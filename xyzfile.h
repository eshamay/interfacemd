#pragma once
#ifndef XYZFILE_H_
#define XYZFILE_H_

#include "vecr.h"
#include "atom.h"
#include <vector>
#include <string>
#include <stdlib.h>


class XYZFile {
	protected:

		Atom_ptr_vec  _atoms;		// The listing of the atoms in the file

		FILE *_file;				// the XYZ file listing all the atom coordinates
		std::string _path;

		int _currentstep, _firstStep, _lastStep, _numSteps,
				_numatoms;				// total number of centers to process from the file for the frame

		bool _initialized;				// To tell wether or not a file has been loaded

		void _FindSteps ();			// assuming that each timestep is headed by an "i = ..." line, then we can load info on the first, last, and total timesteps

	public:

		XYZFile (std::string path);
		~XYZFile ();

		static void WriteXYZ (Atom_ptr_vec& atoms);

		// Various control functions
		// see LoadFirst for the init arg
		void LoadNext ();

		// If the atomlist has been initialized already, then send true, otherwise to form a new list, send false.
		void LoadFirst ();
		void Seek (int step);

		// output functions
		Atom_ptr_vec& Atoms () { return _atoms; }
		Atom_it begin () const { return _atoms.begin(); }
		Atom_it end () const { return _atoms.end(); }

		int Current () const { return _currentstep; }
		int NumSteps () const { return _numSteps; }
		size_t size () const { return _atoms.size(); }

		AtomPtr operator[] (int index) { return _atoms[index]; }

};

#endif
