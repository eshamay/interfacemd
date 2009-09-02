#ifndef _DENSITYTEST_H_
#define _DENSITYTEST_H_

#include "../watersystem.h"

using namespace std;

//#define AVG
//#define RESTART

class DensityAnalyzer : public WaterSystem<AmberSystem> {

protected:

	vector< vector<int> > histo;	// histogram for atom particle numbers arranged by atom names
	vector<string> atomNames;

	void OutputData ();

public:

	vector<int> AtomDensity (string const atomname);
	void SystemDensities ();
	DensityAnalyzer (const int argc, const char **argv, const WaterSystemParams& params);

};

#endif
