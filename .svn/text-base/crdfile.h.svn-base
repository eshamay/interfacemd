#ifndef CRDFILE_H_
#define CRDFILE_H_

#include <stdio.h>
#include <vector>
#include <string>
#include "vecr.h"
#include "utility.h"

class CRDFile {
	
	FILE 			*_file;
	std::vector<VecR>	_coords;	// atomic coordinates
	int 			_size;		// number of atoms in the system
	VecR			_dims;		// Dimensions of the system (box size)
	int 			_frame;		// The current frame (number of timesteps processed)

	bool			_eof;		// end of file marker for the coord file

public:
	
	CRDFile (string crdpath, int size);

	// Various control functions
	void LoadFirst ();
	void LoadNext ();

	// output functions
	const std::vector<VecR>& Coords () const { return _coords; }
	int size () 	const { return _size; }
	const VecR& Dims () const { return _dims; }

	bool eof () 	const { return _eof; }
	int Current () 	const { return _frame; }

	VecR& operator[] (int index) { return _coords[index]; }
};

#endif
