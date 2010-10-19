#include "amber-morita2002.h"

namespace morita {

	void AmberMorita2008Analysis::SelectAnalysisWaters () {

		// slice the waters in the system according to the center-of-mass location
		/*
		// find the system center of mass
		VecR com = Analyzer<AmberSystem>::CenterOfMass (analysis_wats.begin(), analysis_wats.end());
		double cutoff = com[WaterSystem<AmberSystem>::axis] + 20.0;

		std::cout << "cutoff set to: " << cutoff << std::endl;
		*/

		// sort the waters in the system according to position along the main axis
		std::sort(analysis_wats.begin(), analysis_wats.end(), Analyzer<AmberSystem>::molecule_position_pred(Atom::O));

		/*
		for (Morita_it it = analysis_wats.begin(); it != analysis_wats.end(); it++) {
			(*it)->GetAtom(Atom::O)->Print();
		}
		*/


		// only use waters found above a particular location in the slab (i.e. center of mass)
		//std::pair<double,double> slice = std::make_pair(cutoff,WaterSystem<AmberSystem>::posmax);
		//WaterSystem<AmberSystem>::SliceWaters<MoritaH2O_ptr> (analysis_wats, slice);

		//std::cout << analysis_wats.size() << " waters to be analyzed" << std::endl;
		int numAnalysisWaters = WaterSystem<AmberSystem>::SystemParameterLookup ("analysis.morita2002.number-of-analysis-waters");
		analysis_wats.erase (analysis_wats.begin(), analysis_wats.end() - numAnalysisWaters - 1500);
		analysis_wats.erase (analysis_wats.begin() + 50, analysis_wats.end());
		//std::cout << analysis_wats.size() << " waters to be analyzed" << std::endl;

		return;
	}

} // namespace morita





// used to calculate the SFG spectrum based on the morita/hynes 2008 method
int main (int argc, char **argv) {

	Analyzer<AmberSystem> analyzer;
	morita::AmberMorita2008Analysis analysis;
	analyzer.SystemAnalysis(analysis);

	return 0;
}

