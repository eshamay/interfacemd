#pragma once
#ifndef XYZFILE_H_
#define XYZFILE_H_

#include "mdfiles.h"

class XYZFile : public md_files::CoordinateFile {

	void _FindSteps ();			// assuming that each timestep is headed by an "i = ..." line, then we can load info on the first, last, and total timesteps

	public:

	XYZFile (std::string path) : CoordinateFile (path) { }

	// Various control functions
	// see LoadFirst for the init arg
	void LoadNext ();

	// If the atomlist has been initialized already, then send true, otherwise to form a new list, send false.
	void LoadFirst ();

};

#endif
