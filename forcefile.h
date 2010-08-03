#pragma once
#ifndef FORCEFILE_H_
#define FORCEFILE_H_

#include "vecr.h"
#include <iostream>
#include <stdio.h>
#include <string>


class ForceFile {

	FILE 			*_file;
	VecR_vec	 	_forces;	// atomic coordinates
	int 			_size;		// number of atoms in the system
	int 			_frame;		// The current frame (number of timesteps processed)

	bool			_eof;		// end of file marker for the force file
	bool			_loaded;
	bool 			_set;

public:

	ForceFile ();
	ForceFile (std::string const forcepath, int const c_size);
	~ForceFile ();

	// Various control functions
	void LoadFirst ();
	void LoadNext ();

	VecR_it begin () const { return _forces.begin(); }
	VecR_it end () const { return _forces.end(); }

	// output functions
	const VecR_vec& Forces () const { return _forces; }
	int size () 	const { return _size; }

	bool eof () 	const { return _eof; }
	int Current () 	const { return _frame; }
	bool Loaded ()	const { return _loaded; }

};

#endif
