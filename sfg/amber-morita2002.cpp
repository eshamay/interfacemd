#include "amber-morita2002.h"

namespace morita {

	void AmberMorita2008Analysis::SelectAnalysisWaters () {
		// find the system center of mass
		VecR com = Analyzer<AmberSystem>::CenterOfMass (analysis_wats.begin(), analysis_wats.end());
		double cutoff = com[WaterSystem<AmberSystem>::axis] + 5.0;

		// only use waters found above a particular location in the slab (i.e. center of mass)
		std::pair<double,double> slice = std::make_pair(cutoff,WaterSystem<AmberSystem>::posmax);
		WaterSystem<AmberSystem>::SliceWaters<MoritaH2O_ptr> (analysis_wats, slice);

		return;
	}

	void AmberMorita2008Analysis::SetAnalysisWaterDipoleMoments () {
		std::for_each (analysis_wats.begin(), analysis_wats.end(), std::mem_fun(&MoritaH2O::SetDipoleMoment));
		return;
	}

} // namespace morita





// used to calculate the SFG spectrum based on the morita/hynes 2008 method
int main (int argc, char **argv) {

	Analyzer<AmberSystem> analyzer;
	morita::AmberMorita2008Analysis analysis ("amber-morita2002.dat");
	analyzer.SystemAnalysis(analysis);

	return 0;
}

