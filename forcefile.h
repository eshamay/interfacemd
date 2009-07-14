#ifndef FORCEFILE_H_
#define FORCEFILE_H_

#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include "vecr.h"
#include "utility.h"


class ForceFile {

	FILE 			*_file;
	std::vector<VecR>	_forces;	// atomic coordinates
	int 			_size;		// number of atoms in the system
	int 			_frame;		// The current frame (number of timesteps processed)

	bool			_eof;		// end of file marker for the force file
	bool			_loaded;

public:

	ForceFile ();
	ForceFile (std::string const forcepath, int const c_size);
	~ForceFile ();

	// Various control functions
	void LoadFirst ();
	void LoadNext ();

	// output functions
	const std::vector<VecR>& Forces () const { return _forces; }
	int size () 	{ return _size; }

	bool eof () 	{ return _eof; }
	int Current () 	{ return _frame; }
	bool Loaded ()	{ return _loaded; }

	VecR& operator[] (int index) { return _forces[index]; }
};

#endif
