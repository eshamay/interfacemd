#ifndef _DENSITYTEST_H_
#define _DENSITYTEST_H_

#include "../ambersystem.h"
#include "../utility.h"

using namespace std;

//#define AVG
//#define DEBUG

#define PRMTOP	"prmtop"
#define MDCRD	"mdcrd"
#define FORCE	""

#define TIMESTEPS		31800
#define OUTPUT_FREQ		100
#define AXIS			y

#define BINSIZE			0.1

#ifdef AVG
	#define START			-40.0
	#define END				40.0

#else
	#define START			-5.0
	#define END				100.0
#endif

#define PBCFLIP			15.0

#define INT_HIGH	82.479
#define INT_LOW		35.005

class DensityAnalyzer {

private:
	
	FILE *  _output;
	AmberSystem	_sys;

	int 	_step, _steps;		// number of timesteps
	int 	_posbins;		// size of the histogram
	double 	_start;		// starting position for the histogram
	double 	_end;		// ending position
	double 	_binsize;	// size of the bins of the histrogram
	coord 	_axis;		// axis we're going to monitor

	double int_high, int_low, middle;

	vector< vector<int> > _density;
	vector<string> _atomNames;

	void _PrintToFile ();
	void _PrintStatus (int step);

public:

	DensityAnalyzer (char * argv[], int const numAtoms, int const numSteps, double const start, double const end, double const binsize, coord axis);
	vector<int> AtomDensity (string const atomname);
	void SystemDensities ();
	void Debug (string msg) const;

};


#endif
