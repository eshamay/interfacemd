#include "gromacs-morita2002.h"

namespace morita {

	template <class T>
	void GMXMorita2008Analysis<T>::SelectAnalysisWaters () {

		// sort the waters in the system according to position along the main axis
		std::sort(this->analysis_wats.begin(), this->analysis_wats.end(), Analyzer< gromacs::GMXSystem<T> >::molecule_position_pred(Atom::O));

		//std::cout << analysis_wats.size() << " waters to be analyzed" << std::endl;
		int numAnalysisWaters = WaterSystem< gromacs::GMXSystem<T> >::SystemParameterLookup ("analysis.morita2002.number-of-analysis-waters");
		this->analysis_wats.erase (this->analysis_wats.begin(), this->analysis_wats.end() - numAnalysisWaters);
		//this->analysis_wats.erase (this->analysis_wats.begin(), this->analysis_wats.end() - numAnalysisWaters - 200);
		//this->analysis_wats.erase (this->analysis_wats.begin()+numAnalysisWaters, this->analysis_wats.end());
		//std::cout << analysis_wats.size() << " waters to be analyzed" << std::endl;
		return;
	}

	template <class T>
	void GMXMorita2008Analysis<T>::SetAnalysisWaterDipoleMoments () {
		this->MoritaH2ODipoles ();
		return;
	}

	template <class T>
	void GMXMorita2008Analysis<T>::SetAnalysisWaterPolarizability () {
		this->MoritaH2OPolarizabilities();
		return;
	}

} // namespace morita





// used to calculate the SFG spectrum based on the morita/hynes 2008 method
int main (int argc, char **argv) {

	Analyzer<gromacs::GMXSystem<gromacs::XTCFile> > analyzer;
	morita::GMXMorita2008Analysis<gromacs::XTCFile> analysis;
	analyzer.SystemAnalysis(analysis);

	return 0;
}

