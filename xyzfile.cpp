#include "xyzfile.h"

XYZFile::XYZFile (std::string path) {
	_initialized = false;

	_path = path;

	_file = (FILE *)NULL;
	_file = fopen (path.c_str(), "r");
	if (_file == (FILE *)NULL)
	{
		printf ("Couldn't load the XYZ coordinate file file:: %s\n", path.c_str());
		exit(1);
	}

	// find the number of timesteps in the file
	this->_FindSteps();
	// Initialize the atoms
	this->LoadFirst();
}

XYZFile::~XYZFile () {
	fclose (_file);
	RUN(_atoms)
		delete _atoms[i];
}

void XYZFile::LoadNext () {

	char line[1000];

	// first process the frame's header
	fscanf (_file, " %d", &_numatoms);
	//fscanf (_file, " i = %d", &_currentstep);			// some output xyz files have timestep data added
	fgets (line, 1000, _file);	// clear the line out
	fgets (line, 1000, _file);	// clear the line out

	double X, Y, Z;
	char name[10];

	// if we haven't already done so, let's create all the atoms we'll need
	if (!_initialized) {
		_atoms.clear();		// first empty out the list
		for (int i = 0; i < _numatoms; i++) {
			_atoms.push_back (new Atom());
		}
	}

	for (int atom = 0; atom < _numatoms; atom++) {
		// now parse each line's information into the atoms
		fscanf (_file, " %s %lf %lf %lf", name, &X, &Y, &Z);

		// if we haven't already done so, let's create all the atoms we'll need
		if (!_initialized) {
			_atoms[atom]->Name (std::string(name));
			_atoms[atom]->ID (atom);
			_atoms[atom]->SetMass();
			_atoms[atom]->SetCharge();
		}

		// finally we set the position of each atom for the timestep
		_atoms[atom]->Position(X, Y, Z);

	}
	_initialized = true;
	_currentstep++;

return;
}

/*
// Take the next line of the file, and parse the atom information from it, compiling it all into the atom given
Atom XYZFile::_ParseAtom (Atom * pAtom) {

	double x, y, z;

	fscanf (_file, " %s %lf %lf %lf", name, &x, &y, &z);    // parse the info from the line into holders
	//VecR position (x, y, z);
	pAtom->Name (string(name));
	Atom newAtom (name, VecR (x, y, z));
	newAtom.ID(_atoms.size());
	_atoms.push_back (newAtom);

return (newAtom);
}
*/

// grab the number of timesteps in the file
void XYZFile::_FindSteps () {


	// form the command that will get us our number of timesteps using shell commands (should be a lot faster)
	std::string string1 = "cat ";	// Here we'll grab the first timestep
	std::string string2 = " | grep i | head -n 1 | awk '{print $3}'";
	//cout << string1 + _path + string2 << endl;
	FILE * fp = popen((string1 + _path + string2).c_str(), "r");

	fscanf(fp, "%d,", &_firstStep);
	pclose(fp);

	string1 = "tac ";	// now find the last timestep
	fp = popen((string1 + _path + string2).c_str(), "r");
	// and load it up
	fscanf(fp, "%d,", &_lastStep);
	pclose(fp);

	_numSteps = _lastStep - _firstStep;

//	char line[1000];
//	while (!feof(_file)) {
//		fgets (line, 1000, _file);
//		sscanf (line, " i = %d", &_laststep);
//	}

return;
}

void XYZFile::LoadFirst() {
	rewind (_file);
	_currentstep = 0;
	this->LoadNext();
}

void XYZFile::Seek (int step) {

	rewind (_file);

	while (_currentstep != step) {
		this->LoadNext();
	}
}
