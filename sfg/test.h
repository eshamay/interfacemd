#ifndef MORITA_H_
#define MORITA_H_

#include "../utility.h"
#include "../analysis.h"

class Test : public Analyzer<AmberSystem>
{

  public:
	Test (WaterSystemParams& params);

  private:

	void Setup ();
	void Analysis () { return; }
	void DataOutput (const unsigned int timestep) { return; }
	void PostAnalysis () { return; }

};

#endif
