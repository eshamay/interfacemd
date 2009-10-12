#ifndef ANALYSIS_H_
#define ANALYSIS_H_

#include "watersystem.h"

typedef vector<int>(*AxisFunc)(const void *);

class Analyzer : public WaterSystem<AmberSystem> {

  public:
	Analyzer (
		void * analysis_params,
		WaterSystemParams& params
		);
	virtual ~Analyzer ();

	void SystemAnalysis ();

	// create a histogram of the angle between a given molecular axis vector (determined by the axisFunc) and the system's ref_axis. The molecule is chosen by the residue name. The molecules must themselves have the functions for determining the molecular axis vector.
	vector<int> Molecular_Axis_Orientation_Histogram (const string name, VecR (*axisFunc)());

  protected:
	void _OutputHeader () const;
	void _OutputStatus (const int timestep) const;

	// parameters for this particular analysis
	void * _ap;

	// function to perform some initial setup before the main analysis loop
	virtual Setup () const;
	// the main analysis function - this gets run on every timestep
	virtual Analysis () const;
	// something to do after the loop (normalization, etc.) - done after the last timestep
	virtual PostAnalysis () const;
	// A function that defines how data is output from the program
	virtual DataOutput (const int timestep) const;
};

#endif
