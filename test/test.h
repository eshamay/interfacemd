#ifndef TEST_H_
#define TEST_H_

#include "../utility.h"
#include "../analysis.h"

class Tester : public Analyzer<AmberSystem>
{

  public:
	Tester (WaterSystemParams& params);

  private:

	void Setup ();
	void Analysis ();
	void DataOutput (const unsigned int timestep);
	void PostAnalysis () { return; }

};

#endif
