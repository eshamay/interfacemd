#pragma once
#ifndef TOPFILE_H_
#define TOPFILE_H_

#include <cstdlib>
#include <string>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>


class TOPFile {

protected:

	FILE * _topfile;			// The associated topology file output by AMBER

	std::vector<std::string> _atomnames;		// The listing of the atoms in the file
	std::vector<double> _masses;			// atomic masses
	std::vector<double> _charges;		// atomic charges

	std::vector<std::string> _molnames;		// molecule names
	std::vector<int>	_molpointers;	// prmtop pointers to each molecule (the location in the prmtop file, not c-style)
	std::vector<int>	_molsizes;		// the number of atoms in each molecule

	int _numAtoms;
	int _numMols;

public:

	TOPFile (std::string path);
	TOPFile (const TOPFile& topfile);
	TOPFile ();
	~TOPFile ();

	// Various control functions
	void FindFlag (std::string flag);		// To search through the file and find a particular section for parsing
	void LoadSection(std::string flag);

	// output functions
	FILE * File () const { return _topfile; }
	std::vector<std::string>& AtomNames () { return _atomnames; }
	std::vector<double>& Masses () 	{ return _masses; }
	std::vector<double>& Charges () 	{ return _charges; }
	std::vector<std::string>& MolNames () { return _molnames; }
	std::vector<int>& MolPointers () { return _molpointers; }
	std::vector<int>& MolSizes ()	{ return _molsizes; }

	int NumAtoms () const { return _numAtoms; }
	int NumMols () const { return _numMols; }

};

#endif
