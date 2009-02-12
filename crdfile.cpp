#include "crdfile.h"

CRDFile::CRDFile (std::string crdpath, int size) : _size(size) {

	// first load up the file given the path
	_file = (FILE *)NULL;
	_file = fopen64 (crdpath.c_str(), "r");
	if (_file == (FILE *)NULL) {
		printf ("Couldn't load the crdfile %s\n", crdpath.c_str());
	}

	_eof = false;

	char str[1000];	
	fgets(str, 1000, _file);		// first frame's header
	this->LoadFirst ();	// load the first frame of the file

}

void CRDFile::LoadNext () {

	this->_coords.clear();
	double x, y, z;
	//char line[1000];

	// grab each coordinate vector for each atom until the size of the system (# of atoms) is processed
	for (int i = 0; i < _size; i++) {
		if (fscanf (_file, " %lf %lf %lf", &x, &y, &z) == EOF) {
			// if we've reached the end of the file, then let it be known
			_eof = true;
		}
		else {
			_coords.push_back(VecR (x, y, z));
		}
	}

	// process the next frame's header line (grab the box dimensions)
	fscanf (_file, " %lf %lf %lf", &x, &y, &z);
	_dims = VecR (x, y, z);
	
	_frame++;

return;
}

void CRDFile::LoadFirst() {
	rewind (_file);
	this->LoadNext();
	_frame = 1;
}
