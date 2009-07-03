#include "crdfile.h"

CRDFile::CRDFile (std::string const crdpath, int const c_size) :
	_file((FILE *)NULL),
	_size(c_size),
	_dims(),
	_frame(0),
	_eof(true)
{

	// first load up the file given the path
	_file = fopen64 (crdpath.c_str(), "r");
	if (_file == (FILE *)NULL) {
		printf ("Couldn't load the crdfile %s\n", crdpath.c_str());
		exit(1);
	}

	_eof = true;

	char str[1000];
	fgets(str, 1000, _file);		// first frame's header
	this->LoadFirst ();	// load the first frame of the file

}

CRDFile::~CRDFile () {
	fclose(_file);
}

void CRDFile::LoadNext () {

	_coords.resize(_size, VecR());
	double x, y, z;
	//char line[1000];

	// grab each coordinate vector for each atom until the size of the system (# of atoms) is processed
	for (int i = 0; i < _size; i++) {
		if (fscanf (_file, " %lf %lf %lf", &x, &y, &z) == EOF) {
			// if we've reached the end of the file, then let it be known
			_eof = true;
		}
		else {
			_coords[i].Set(x,y,z);
		}
	}

	// process the next frame's header line (grab the box dimensions)
	fscanf (_file, " %lf %lf %lf", &x, &y, &z);
	_dims.Set(x,y,z);

	_frame++;

return;
}

void CRDFile::LoadFirst() {
	rewind (_file);
	this->LoadNext();
	_frame = 1;
}
