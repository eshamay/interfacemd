#include "forcefile.h"

ForceFile::ForceFile (std::string const forcepath, int const c_size) :
	_file((FILE *)NULL), _size(c_size), _eof(true), _loaded(false), _set(false)
{

	// first load up the file given the path
	_file = fopen64 (forcepath.c_str(), "r");
	if (_file != (FILE *)NULL) {
		_loaded = true;

		_eof = false;

		char str[1000];
		fgets(str, 1000, _file);		// first frame's header
		this->LoadFirst ();	// load the first frame of the file
	}
}

	ForceFile::~ForceFile () {
		if (_loaded)
			fclose(_file);
	}

void ForceFile::LoadNext () {

	if (!_set) {
		_forces.resize(_size, VecR());
		_set = true;
	}
	double x, y, z;
	//	char line[1000];

	// grab each force vector for each atom until the size of the system (# of atoms) is processed
	for (VecR_it_non_const it = _forces.begin(); it != _forces.end(); it++) {
		if (fscanf (_file, " %lf %lf %lf", &x, &y, &z) == EOF) {
			// if we've reached the end of the file, then let it be known
			_eof = true;
		}
		else {
			it->Set(x,y,z);
		}
	}

	++_frame;

	return;
}

void ForceFile::LoadFirst() {
	rewind (_file);
	_set = false;
	this->LoadNext();
}
