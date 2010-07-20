#ifndef BOXFILLER_H_
#define BOXFILLER_H_

#include "pdbfile.h"
#include "vecr.h"
#include <vector>
#include <iostream>
#include <math.h>

#define		PARAMS		string(argv[1])
#define		PDBFILE		string(argv[2])
#define 	SPACING		atoi(argv[3])
#define 	RANDOM  	double(rand())/RAND_MAX

using namespace std;

class BoxFiller {

private:
	PDBFile _pdb;		// pdb file containing all the atoms/molecules
	FILE * _params;
	double	_spacing;

	VecR	_boxSize;	// system size

	vector<string>	_residueNames;
	vector<int>		_residueNum;
	vector<Molecule>	_mols;

	void _InitParams ();

	void _FillBoxRandom ();
	void _FillBoxLattice ();

public:

	BoxFiller (string paramfile, string pdbfile, double spacing);
	~BoxFiller ();

};

#endif
