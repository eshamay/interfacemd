#ifndef XYZSFG_H_
#define XYZSFG_H_

#include "../analysis.h"
#include "fftw3.h"

class XYZSFGAnalyzer : public Analyzer<XYZSystem>
{

  public:
	XYZSFGAnalyzer (WaterSystemParams& params);

  private:

	void Setup ();
	void Analysis ();
	void DataOutput (const unsigned int timestep);
	void PostAnalysis () { return; }

	VecR_vec	_M;

};

#endif

