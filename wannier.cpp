#include "wannier.h"

WannierFile::WannierFile ()
  :
	_file((FILE *)NULL),
	_loaded (false)
{ }


WannierFile::WannierFile (std::string wannierpath) 
  :
	_file((FILE *)NULL)
{
  // first load up the file given the path
  _file = (FILE *) NULL;
  _file = fopen64 (wannierpath.c_str(), "r");
  if (!_file) {
	printf ("WannierFile c-tor - Couldn't open the wannier file specified:\n\t%s\n", wannierpath.c_str());
	exit(1);
  }
  else {
	_loaded = true;
	_eof = false;
	this->LoadFirst ();	// load the first frame of the file
  }

}

WannierFile::~WannierFile () {
  if (_loaded) fclose(_file);
  return;
}

void WannierFile::LoadNext () {

  if (_loaded) {

	_coords.clear();

	// grab info from the header: number of centers, and frame number
	fscanf (_file, " %d 1_%d ", &_size, &_frame);

	double a, b, c;

	// grab each coordinate vector for each wannier center until the size of the system is processed
	for (int i = 0; i < _size; i++) {
	  //printf ("% 8.4f % 8.4f % 8.4f\n", a,b,c);
	  //if (fscanf (_file, " %*s %lf %lf %lf %*lf ", &a, &b, &c) == EOF) _eof=true;
	  if (fscanf (_file, " %*s %lf %lf %lf %*f %*f %*f ", &a, &b, &c) == EOF) _eof=true;
	  _coords.push_back(VecR (a, b, c));
	}

  }

  return;
}

void WannierFile::LoadFirst() {
  if (_loaded) {
	this->LoadNext();
  }
}
