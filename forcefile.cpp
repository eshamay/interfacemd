#include "forcefile.h"

ForceFile::ForceFile (string forcepath, int size) : _size(size) {
	_loaded = false;
	
	// first load up the file given the path
	_file = (FILE *)NULL;
	_file = fopen64 (forcepath.c_str(), "r");
	if (_file == (FILE *)NULL) {
		//printf ("Couldn't load the force file %s\n", forcepath.c_str());
		_loaded = false;
	}
	else {
		_loaded = true;

		_eof = false;

		char str[1000];	
		fgets(str, 1000, _file);		// first frame's header
		this->LoadFirst ();	// load the first frame of the file
	}
}

void ForceFile::LoadNext () {

	_forces.clear();
	double x, y, z;
//	char line[1000];

	// grab each force vector for each atom until the size of the system (# of atoms) is processed
	for (int i = 0; i < _size; i++) {
		if (fscanf (_file, " %lf %lf %lf", &x, &y, &z) == EOF) {
			// if we've reached the end of the file, then let it be known
			_eof = true;
		}
		else {
			_forces.push_back(VecR (x, y, z));
		}
	}

	_frame++;

return;
}

void ForceFile::LoadFirst() {
	rewind (_file);
	this->LoadNext();
	_frame = 1;
}
