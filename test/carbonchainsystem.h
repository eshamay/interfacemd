#ifndef	CARBONCHAINSYSTEM_H_
#define	CARBONCHAINSYSTEM_H_

#include "../carbonchain.h"
#include "../decane.h"
#include "../analysis.h"
#include "../utility.h"

typedef struct {
	string mol_name;
	vector<int> histogram; 					// total running histogram
} AnalysisParams;

typedef Decane molecule_t;

class CarbonChainSystem : public Analyzer<AnalysisParams> {

	public:

		CarbonChainSystem 
			(
			 AnalysisParams& analysis_params,
			 WaterSystemParams& params
			);

		void Setup ();
		void Analysis ();
		void DataOutput (const int timestep);
		void PostAnalysis ();

};

#endif
