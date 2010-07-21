#include "test.h"


Tester::Tester (WaterSystemParams& wsp)
  :	
	Analyzer<XYZSystem> (wsp)
{
  return;
}

void Tester::Setup () {

  this->sys->SetReparseLimit(1);

  return;
}

void Tester::Analysis () {

  LoadAll();

  Mol_it so2 = md_utility::FindByName(sys_mols.begin(), sys_mols.end(), "so2");
  //printf ("the so2 is:\n");
  //(*so2)->Print();

  bondgraph::distance_pair dt = this->sys->graph.ClosestAtom(*so2);
  min_distances.push_back(dt);


  return;

}

void Tester::DataOutput () {
  rewind(output);

  for (bondgraph::distance_vec::const_iterator it = min_distances.begin(); it != min_distances.end(); it++) {
	fprintf (output, "% 8.3f % 8d %8s/%s\n", it->first, it->second->MolID(), it->second->Name().c_str(), it->second->Residue().c_str());
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

