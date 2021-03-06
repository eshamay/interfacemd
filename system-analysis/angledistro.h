#ifndef ANGLEDISTRO_H_
#define ANGLEDISTRO_H_

#include <stdio.h>
#include <stdlib.h>

#include "../utility.h"
#include "../ambersystem.h"

#define PRMTOP  "prmtop"
#define MDCRD   "mdcrd"
#define FORCE   "mdvel"

#define AXIS	y

#define OUTPUT_FILE	"angledistro.dat"

#define	ZMAX	100.0
#define	ZMIN	-5.0
#define	ZRES	0.1

#define AMAX	1.0
#define	AMIN	-1.0
#define ARES	0.02

class AngleDistro {

private:
	AmberSystem _sys;
	coord _axis;

	double _zmin, _zmax, _zres;
	int _zbins;

	double _anglemin, _anglemax, _angleres;
	int _anglebins;

	FILE * _output;

	vector< vector<int> > _histo;
	void _ClearHisto ();

	int _timestep;

	void _PrintStatus () const;
	void _PrintOutput () const;
	vector< vector<int> > _WaterStepHistogram ();

public:

	AngleDistro ();
	~AngleDistro ();

	vector< vector<int> >& WaterHistogram (int const timesteps);

};



#endif
