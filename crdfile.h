#ifndef CRDFILE_H_
#define CRDFILE_H_

#include "vecr.h"
#include <stdio.h>
#include <string>

class CRDFile {

	FILE 			*_file;
	VecR_vec		_coords;	// atomic coordinates
	int 			_size;		// number of atoms in the system
	VecR			_dims;		// Dimensions of the system (box size)
	int 			_frame;		// The current frame (number of timesteps processed)

	bool			_eof;		// end of file marker for the coord file
	bool			_set;
	bool			_periodic;	// are periodic boundaries being used

	public:

	CRDFile () { }
	CRDFile (std::string const crdpath, int const c_size, const bool periodic);
	~CRDFile ();

	// Various control functions
	void LoadFirst ();
	void LoadNext ();

	VecR_it begin () const { return _coords.begin(); }
	VecR_it end () const { return _coords.end(); }
	// output functions
	const VecR_vec& Coords () const { return _coords; }
	int size () 	const { return _size; }
	const VecR& Dims () const { return _dims; }

	bool eof () 	const { return _eof; }
	int Current () 	const { return _frame; }

};

#endif
