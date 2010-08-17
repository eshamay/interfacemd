#include "cp2k-morita2002.h"

namespace morita {

  CP2KSFGAnalyzer::CP2KSFGAnalyzer (WaterSystemParams& wsp) : SFGAnalyzer<XYZSystem> (wsp) { return; }


// for the small cp2k systems, just use all the waters
  void CP2KSFGAnalyzer::SelectAnalysisWaters () {
	return;
  }

  void CP2KSFGAnalyzer::SetAnalysisWaterDipoleMoments () {
	std::for_each (analysis_wats.begin(), analysis_wats.end(), MDSystem::CalcDipole);
	return;
  }

} // namespace morita




// used to calculate the SFG spectrum based on the morita/hynes 2002 method
int main (int argc, char **argv) {

  libconfig::Config cfg;
  cfg.readFile("system.cfg");

  std::string filename ("cp2k-morita2002.dat");
  libconfig::Setting &analysis = cfg.lookup("analysis");
  analysis.add("filename", libconfig::Setting::TypeString) = filename;

  WaterSystemParams wsp (cfg);

  morita::CP2KSFGAnalyzer sfg (wsp);

  sfg.SystemAnalysis ();

  return 0;
}
