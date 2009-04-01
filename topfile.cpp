#include "topfile.h"

TOPFile::TOPFile (std::string path) :
	_topfile((FILE *)NULL)
	{

	// Before anything, let's load the file!
	_topfile = fopen(path.c_str(), "r");
	if (_topfile == (FILE *)NULL) {
		std::cout << "Error opening the topfile " << path << std::endl;
		return;
	}

	// First thing is to go through the topology file and parse out all the atom information associated with it.

	// the "POINTERS" section stores various system-wide values as desccribed in the topology file specification in AMBER docs
	this->FindFlag("POINTERS");
	fscanf (_topfile, " %d", &_numAtoms);	// Here we grab the number of atoms in the system

	this->LoadSection("ATOM_NAME");
	this->LoadSection("MASS");
	this->LoadSection("CHARGE");
	this->LoadSection("RESIDUE_LABEL");
	this->LoadSection("RESIDUE_POINTER");
	this->LoadSection("ATOMS_PER_MOLECULE");
	
	_numMols = _molnames.size();
return;
}

TOPFile::TOPFile (const TOPFile& topfile) :
	_topfile(topfile._topfile),
	_atomnames(topfile._atomnames),
	_masses(topfile._masses),
	_charges(topfile._charges),
	_molnames(topfile._molnames),
	_molpointers(topfile._molpointers),
	_molsizes(topfile._molsizes),
	_numAtoms(topfile._numAtoms),
	_numMols(topfile._numMols)
{}

TOPFile::TOPFile () {
}

TOPFile::~TOPFile () {
	fclose(_topfile);
}

void TOPFile::FindFlag (std::string flag) {
	rewind (_topfile);
	char str[1000] = "";

	// run through the file until the flag is found
	while (strcmp(str, flag.c_str())) {
		fgets (str, 1000, _topfile); 		// load the string up 
		sscanf (str, " %*s %s", str);
	}
	// then strip the next line which holds a 'format' for the section
	fgets (str, 1000, _topfile); 		// load the string up 

return;
}

/* 
The topfile uses the %FLAG ATOM_NAME to mark the beginning of the atom-names section. After that point the atom names are all listed sequentially. A very easy parsing job.
*/
// As with ATOM_NAME, the rest of the sections will parse out various properties of each atom or molecule
void TOPFile::LoadSection(std::string flag) {
	this->FindFlag (flag);

	char value[1000] = "";

	// run through each line of the section and pull out the values of the atoms
	while (strcmp(value, "%FLAG")) {
		fscanf (_topfile, " %s", value);

		if (strcmp(value, "%FLAG")) {
			if (flag == "ATOM_NAME") {
				_atomnames.push_back(value);
			}
			if (flag == "MASS") {
				_masses.push_back(atof(value));
			}
			if (flag == "CHARGE") {
				_charges.push_back(atof(value));
			}
			if (flag == "RESIDUE_LABEL") {
				_molnames.push_back(value);
			}
			if (flag == "RESIDUE_POINTER") {
				_molpointers.push_back(atoi(value));
			}
			if (flag == "ATOMS_PER_MOLECULE") {
				_molsizes.push_back(atoi(value));
			}
		}
	}

return;
}

