#ifndef MORITA_H_
#define MORITA_H_

#include "../utility.h"
#include "../analysis.h"
#include "watersfg.h"

typedef std::vector< complex<double> > Complex_vec;

class SFGAnalyzer : public Analyzer 
{

  public:
	SFGAnalyzer (WaterSystemParams& params);

  private:

	SFGCalculator	sfg;
	Complex_vec Molecular_Beta;			// chi for a given molecule
	Complex_vec TimestepChi;	// chi spectrum of several different molecules collected over an entire timestep
	Complex_vec TotalChi;		// Running total of all the data for several timesteps

	unsigned long numMolsProcessed;
	bool firstmol;
	bool firsttimestep;

	void Setup ();
	void Analysis ();
	void DataOutput (const unsigned int timestep);
	void PostAnalysis ();

	void CollectChi (Complex_vec& newchi, Complex_vec& totalchi);

	void FlipWaters (const coord axis = y);

};

#endif
