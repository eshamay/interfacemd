#include "amber-morita2002.h"

namespace morita {

  AmberSFGAnalyzer::AmberSFGAnalyzer (WaterSystemParams& wsp)
	:
	  SFGAnalyzer<AmberSystem> (wsp)
  { return; }


  void AmberSFGAnalyzer::SelectAnalysisWaters () {
	// find the system center of mass
	VecR com (this->CenterOfMass<MoritaH2O_ptr>(all_wats));
	double cutoff = com[WaterSystem<U>::axis];

	// only use waters found above a particular location in the slab (i.e. center of mass)
	std::pair<double,double> slice = std::make_pair(cutoff,WaterSystem<U>::posmax);
	this->SliceWaters<MoritaH2O_ptr> (analysis_wats, slice);

	return;
  }

  void AmberSFGAnalyzer::SetAnalysisWaterDipoleMoments () {
	std::for_each (analysis_wats.begin(), analysis_wats.end(), std::mem_fun(&MoritaH2O::SetDipoleMoment));
	return;
  }

} // namespace morita




// used to calculate the SFG spectrum based on the morita/hynes 2002 method
int main (int argc, char **argv) {

  libconfig::Config cfg;
  cfg.readFile("system.cfg");

  std::string filename ("amber-morita2002.dat");
  libconfig::Setting &analysis = cfg.lookup("analysis");
  analysis.add("filename", libconfig::Setting::TypeString) = filename;

  WaterSystemParams wsp (cfg);

  morita::AmberSFGAnalyzer sfg (wsp);

  sfg.SystemAnalysis ();

  return 0;
}
