#ifndef TEST_H_
#define TEST_H_

#include "../analysis.h"

class Tester : public Analyzer<XYZSystem>
{

  public:
	Tester (WaterSystemParams& params);

  private:

	void Setup ();
	void Analysis ();
	void DataOutput (const unsigned int timestep);
	void PostAnalysis () { return; }

	std::vector< std::pair<int,double> > min_distances;
	int timestep;
};

#endif

