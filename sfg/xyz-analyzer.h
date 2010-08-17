#ifndef XYZ_H_
#define XYZ_H_

#include "../analysis.h"

class XYZAnalyzer : public Analyzer<XYZSystem>
{

  public:
	XYZAnalyzer (WaterSystemParams& wsp) : Analyzer<XYZSystem> (wsp), timezero(false) { return; }

  private:

	void Setup () { 
	  this->sys->SetReparseLimit(1); 
	  return; 
	}	// Setup

	void Analysis ();
	void DataOutput () { return; }
	void PostAnalysis () { return; }

	bool timezero;
	VecR M_0;
};

#endif

