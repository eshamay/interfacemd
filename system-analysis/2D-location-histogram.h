#ifndef 2D_LOCATION_HISTOGRAM_H_
#define 2D_LOCATION_HISTOGRAM_H_

#include "../utility.h"
#include "../analysis.h"

class InplaneHistogram : public Analyzer {

  public:
	Angle_Position_Analyzer (WaterSystemParams& wsp) :
	  Analyzer(wsp) { }

	// the utility functor for getting all the data accumulated
	APBinner<Water> binner;

	void Setup () {
	  // For now let's look at waters only
	  FindAll ();
	  KEEP_BY_NAME(int_mols, Molecule *, "odn");
	  return;
	}

	void Analysis () {
	  // bin each water's position and angle
	  FOR_EACH(int_mols, this->binner);
	  return;
	}

	// nothing at the moment
	void PostAnalysis () { return; }

	void DataOutput (const unsigned int timestep) {

	  if (!(timestep % (output_freq * 10)))
		this->binner.APBinnerOutput(output);
	  return;
	}

};
#endif
