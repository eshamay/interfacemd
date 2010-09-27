#include "cp2k-morita2002.h"

namespace morita {

	// for the small cp2k systems, just use all the waters
	void CP2KMorita2002Analysis::SelectAnalysisWaters () {
		return;
	}

	void CP2KMorita2002Analysis::SetAnalysisWaterDipoleMoments () {
		this->CalcWannierDipoles ();
		return;
	}

	void CP2KMorita2002Analysis::SetAnalysisWaterPolarizability () {
		this->MoritaH2OPolarizabilities();
		return;
	}

} // namespace morita




// used to calculate the SFG spectrum based on the morita/hynes 2002 method
int main (int argc, char **argv) {

	Analyzer<XYZSystem> analyzer;
	morita::CP2KMorita2002Analysis analysis ("cp2k-morita2002.dat");
	analyzer.SystemAnalysis(analysis);

	return 0;
}
