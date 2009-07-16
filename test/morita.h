#ifndef MORITA_H_
#define MORITA_H_

#include <stdio.h>
#include <stdlib.h>
#include "../watersystem.h"
#include "../watersfg.h"

// if we're running this on an MPI system to do some parallel work:
//#define MPI_SYS
#ifdef	MPI_SYS
	#include "../mpi/mpisys.h"
#endif

typedef std::vector< complex<double> > Complex_vec;

class SFGAnalyzer : public WaterSystem<AmberSystem> {

public:
	SFGAnalyzer (const int argc, const char **argv, const WaterSystemParams& params);
	void Analyze ();

private:

	SFGCalculator	sfg;
	Complex_vec MolChi;			// chi for a given molecule
	Complex_vec TimestepChi;	// chi spectrum of several different molecules collected over an entire timestep
	Complex_vec TotalChi;		// Running total of all the data for several timesteps

	unsigned long numMolsProcessed;

	void OutputData ();
	void CollectChi (Complex_vec& newchi, Complex_vec& totalchi);

};

#endif
