#ifndef MORITA_H_
#define MORITA_H_

#include "../utility.h"
#include "../analysis.h"
#include "watersfg.h"

typedef std::vector< complex<double> > Complex_vec;

class SFGAnalyzer : public Analyzer<AmberSystem>
{

  public:
	SFGAnalyzer (WaterSystemParams& params);

  private:

	static SFGCalculator	sfg;
	static Complex_vec TimestepChi;		// chi spectrum of several different molecules collected over an entire timestep
	static Complex_vec TotalChi;		// Running total of all the data for several timesteps

	static unsigned long numMolsProcessed;
	static bool firstMol;
	static bool firstTimeStep;

	void Setup ();
	void Analysis ();
	void DataOutput (const unsigned int timestep);
	void PostAnalysis ();

	template <typename Iter1, typename Iter2>
	  // Build up running totals of the chi spectra
	static void CollectChi (Iter1 pBegin, Iter1 pEnd, Iter2 pBegin2)
	{
	  typedef typename std::iterator_traits<Iter1>::value_type value_type;

	  std::transform(pBegin, pEnd, pBegin2, pBegin, std::plus<value_type>());
	}

	//static void CollectChi (Complex_vec& newchi, Complex_vec& totalchi);

	void FlipWaters (const coord axis = y);

	class SFGProcessor : public std::unary_function<Molecule *,bool> {
	  private:
		Water * water;
		Complex_vec Molecular_Beta;			// chi for a given molecule

	  public:

		bool operator() (Molecule * mol) {
		  water = static_cast<Water *>(mol);

		  Molecular_Beta.clear();

		  // and then calculate the chi spectrum for the molecule 
		  // 0,2,1 = SSP. S = X and Z axes, P = Y axis
		  Molecular_Beta = SFGAnalyzer::sfg.Beta (*water, 0,2,1);

		  SFGAnalyzer::numMolsProcessed++;

		  // when starting a new timestep...
		  if (SFGAnalyzer::firstMol) {
			SFGAnalyzer::TimestepChi.resize (Molecular_Beta.size(), complex<double> (0.0, 0.0));
			SFGAnalyzer::firstMol = false;
		  }

		  // for the very first timestep...
		  if (SFGAnalyzer::firstTimeStep) {
			SFGAnalyzer::TotalChi.clear();
			SFGAnalyzer::TotalChi.resize (Molecular_Beta.size(), complex<double> (0.0, 0.0));
			SFGAnalyzer::firstTimeStep = false;
		  }

		  // perform the summation for averaging spectra over all the water molecules
		  SFGAnalyzer::CollectChi (SFGAnalyzer::TimestepChi.begin(), SFGAnalyzer::TimestepChi.end(), Molecular_Beta.begin());

		  return true;
		}
	};

	SFGProcessor SFGProcess;
};

#endif
