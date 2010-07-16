#include "test.h"


Tester::Tester (WaterSystemParams& wsp)
  :	
	Analyzer<XYZSystem> (wsp)
{
  return;
}

void Tester::Setup () {

  this->sys->SetReparseLimit(20);

  return;
}

void Tester::Analysis () {

  LoadAll();

  std::vector< mol_distance > distances;
 
  std::string name ("so2");
  Mol_it so2 = md_utility::FindByName(sys_mols.begin(), sys_mols.end(), name);
  //Mol_it so2 = std::find_if (sys_mols.begin(), sys_mols.end(), std::bind2nd(name_pred<MolPtr>(), "so2"));

  for (Mol_it mol = sys_mols.begin(); mol != sys_mols.end(); mol++) {
	if (so2 == mol) continue;
	distances.push_back(std::make_pair(this->sys->Distance(*so2,*mol), *mol));
  }

  md_utility::pair_sort_first (distances.begin(), distances.end());
  min_distances.push_back(distances[0]);

  return;

}

void Tester::DataOutput (const unsigned int timestep) { 

  rewind(output);

  for (int i = 0; i < min_distances.size(); i++) {
	fprintf (output, "% 8.3f % 8d %s\n", min_distances[i].first, min_distances[i].second->MolID(), min_distances[i].second->Name().c_str());
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

