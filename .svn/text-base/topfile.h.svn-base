#ifndef TOPFILE_H_
#define TOPFILE_H_

#include <string>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "utility.h"


class TOPFile {
	
	FILE * _topfile;			// The associated topology file output by AMBER

	std::vector<string> _atomnames;		// The listing of the atoms in the file
	std::vector<double> _masses;			// atomic masses
	std::vector<double> _charges;		// atomic charges

	std::vector<string> _molnames;		// molecule names
	std::vector<int>	   _molpointers;	// prmtop pointers to each molecule (the location in the prmtop file, not c-style)
	std::vector<int>	   _molsizes;		// the number of atoms in each molecule

	int _numAtoms;
	int _numMols;

public:
	
	TOPFile (string path);
	TOPFile (const TOPFile& topfile);
	TOPFile ();
	//~TOPFile ();

	// Various control functions
	void FindFlag (string flag);		// To search through the file and find a particular section for parsing
	void LoadSection(string flag);

	// output functions
	FILE * File () const { return _topfile; }
	std::vector<string> AtomNames () const { return _atomnames; }
	std::vector<double> Masses () 	const { return _masses; }
	std::vector<double> Charges () 	const { return _charges; }
	std::vector<string> MolNames () const { return _molnames; }
	std::vector<int> MolPointers () const { return _molpointers; }
	std::vector<int> MolSizes ()	const { return _molsizes; }

	int NumAtoms () const { return _numAtoms; }
	int NumMols () const { return _numMols; }

};

#endif
