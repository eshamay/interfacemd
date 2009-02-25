#ifndef _DENSITYTEST_H_
#define _DENSITYTEST_H_

#include "../watersystem.h"

using namespace std;

#define AVG
//#define DEBUG

#define PRMTOP	"prmtop"
#define MDCRD	"mdcrd"
#define FORCE	""

#define TIMESTEPS		200000
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

#define INT_HIGH	82.4596
#define INT_LOW		35.98845

class DensityAnalyzer : public WaterSystem {

private:
	
	vector< vector<int> > _density;
	vector<string> _atomNames;

public:

	vector<int> AtomDensity (string const atomname);
	void SystemDensities ();

};

#endif
