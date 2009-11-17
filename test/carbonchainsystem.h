#ifndef	CARBONCHAINSYSTEM_H_
#define	CARBONCHAINSYSTEM_H_

#include "../carbonchain.h"
#include "../decane.h"
#include "../utility.h"
#include "../analysis.h"

typedef VecR (CarbonChain::*Axis_Vector_Fn) ();	// a pointer to an axis vector function to find carbon-chain axes
typedef void (CarbonChain::*FnPtr) ();
typedef void (CarbonChain::*OutputPtr) (const int timestep);

typedef struct {
	string 			mol_name;
	string			atom_name;
	vector<double> 	histogram;
	vector<int>		int_histogram;
	vector< vector<double> >	histogram_2;
	vector< vector<int> >		density_histo_2;
	// A function that returns the long-axis tilt vector for a given molecule
	Axis_Vector_Fn	axisFn;
} AnalysisParams;

typedef Decane mol_t;

class CarbonChainSystem : public Analyzer<CarbonChainSystem, AnalysisParams> {

	public:

		CarbonChainSystem 
			(
			 AnalysisParams& analysis_params,
			 WaterSystemParams& params
			);


		void EmptyFn ();
		void ClearHistogram_FindMols ();
		void MolecularAxisDirection ();
		void AngleHistogramOutput (const int timestep);
		void Histogram_2_Output (const int timestep);
		void Histogram_2_Output_2 (const int timestep);
		void HistogramAveraging ();
		void HistogramAveraging_2 ();
		void MolecularPlane ();
};

#endif
