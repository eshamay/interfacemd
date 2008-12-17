#include "topfile.h"

TOPFile::TOPFile (string path) {

	// Before anything, let's load the file!
	_topfile = (FILE *)NULL;
	_topfile = fopen(path.c_str(), "r");
	if (_topfile == (FILE *)NULL) {
		cout << "Error opening the topfile " << path << endl;
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

TOPFile::TOPFile (const TOPFile& topfile) {
	_topfile = topfile.File();	
	_atomnames = topfile.AtomNames();
	_masses = topfile.Masses();
	_charges = topfile.Charges();
	_molnames = topfile.MolNames();
	_molpointers = topfile.MolPointers();
	_molsizes = topfile.MolSizes();
	_numAtoms = topfile.NumAtoms();
	_numMols = topfile.NumMols();

}

TOPFile::TOPFile () {
}

void TOPFile::FindFlag (string flag) {
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
void TOPFile::LoadSection(string flag) {
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

