#pragma once
#ifndef WANNIER_H_
#define WANNIER_H_

#include "vecr.h"
#include <stdio.h>
#include <string>

#define WANNIER_BOND	0.7


class WannierFile {

protected:

	FILE *			_file;
	VecR_vec		_coords;	// atomic coordinates
	int 			_size;		// number of lines to process from the wannier file per frame
	int 			_frame;		// The current frame (number of timesteps processed)
	int				_ID;

	bool			_eof;		// end of file marker for the coord file
	bool			_loaded;	// tells wether or not the wannier file is being used or not (or if it exists)

public:

	WannierFile (std::string wannierpath);
	~WannierFile ();

	// Various control functions
	void LoadFirst ();
	void LoadNext ();

	// output functions
	const VecR_vec& Coords () const { return _coords; }
	VecR_it begin () const { return _coords.begin(); }
	VecR_it end () const { return _coords.end(); }

	unsigned int size () 	const { return _size; }

	bool eof () 	const { return _eof; }		// have we reached the end of the file?
	bool Loaded ()	const { return _loaded; }	// find out if the file is loaded/exists
	int Current () 	const { return _frame; }

};

#endif
