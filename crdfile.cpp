#include "crdfile.h"

CRDFile::CRDFile (std::string const crdpath, int const c_size) :
	_file((FILE *)NULL), _size(c_size), _dims(), _frame(0), _eof(true), _set(false)
{

	// first load up the file given the path
	_file = fopen64 (crdpath.c_str(), "r");
	if (_file == (FILE *)NULL) {
		printf ("Couldn't load the crdfile %s\n", crdpath.c_str());
		exit(1);
	}

	_eof = false;

	char str[1000];
	fgets(str, 1000, _file);		// first frame's header
	this->LoadFirst ();	// load the first frame of the file

}

CRDFile::~CRDFile () {
	fclose(_file);
}

void CRDFile::LoadNext () {

	if (!_set) {
		_coords.resize(_size, VecR());
		_set = true;
	}
	double x, y, z;

	// grab each coordinate vector for each atom until the size of the system (# of atoms) is processed
	for (VecR_it_non_const it = _coords.begin(); it != _coords.end(); it++) {
		if (fscanf (_file, " %lf %lf %lf", &x, &y, &z) == EOF) {
			// if we've reached the end of the file, then let it be known
			_eof = true;
		}
		else {
			it->Set(x,y,z);
		}
	}

	// process the next frame's header line (grab the box dimensions)
	fscanf (_file, " %lf %lf %lf", &x, &y, &z);
	_dims.Set(x,y,z);

	++_frame;

	return;
}


void CRDFile::LoadFirst() {
	rewind (_file);
	_set = false;
	this->LoadNext();
}
