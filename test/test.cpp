#include "test.h"


Tester::Tester (WaterSystemParams& wsp)
  :	
	Analyzer<XYZSystem> (wsp), timestep(0)
{
  return;
}

void Tester::Setup () {

  this->sys->SetReparseLimit(100);

  return;
}

void Tester::Analysis () {

  LoadAll();

  std::vector<double> distances;
  Mol_it so2 = std::find_if (sys_mols.begin(), sys_mols.end(), std::bind2nd(name_pred<MolPtr>(), "so2"));
  //(*so2)->Print();

  /*
  for (Mol_it mol = sys_mols.begin(); mol != sys_mols.end(); mol++) {
	if (so2 == mol) continue;
	distances.push_back(this->sys->Distance(*so2,*mol));
  }

  std::sort(distances.begin(), distances.end());
  min_distances.push_back(std::make_pair(timestep++,distances[0]));
  */

  return;

}

void Tester::DataOutput (const unsigned int timestep) { 

  rewind(output);

  for (int i = 0; i < min_distances.size(); i++) {
	fprintf (output, "% 8d % 8.3f\n", min_distances[i].first, min_distances[i].second);
  }

  return; 
}



int main () {

  libconfig::Config cfg;
  cfg.readFile("system.cfg");

  std::string filename = cfg.lookup("analysis.test.filename");
  libconfig::Setting &analysis = cfg.lookup("analysis");
  analysis.add("filename", libconfig::Setting::TypeString) = filename;

  WaterSystemParams wsp (cfg);

  Tester test (wsp);

  test.SystemAnalysis ();

  return 0;
}

