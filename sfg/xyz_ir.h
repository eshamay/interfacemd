#ifndef XYZSFG_H_
#define XYZSFG_H_

#include "../analysis.h"

class XYZSFGAnalyzer : public Analyzer<XYZSystem>
{

  public:
	XYZSFGAnalyzer (WaterSystemParams& params);

  private:

	void Setup ();
	void Analysis ();
	void DataOutput ();
	void PostAnalysis () { return; }

	VecR_vec	_M;

};

#endif

