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

  Mol_it so2 = md_utility::FindByName(sys_mols.begin(), sys_mols.end(), "so2");
  printf ("the so2 is:\n");
  (*so2)->Print();

  bondgraph::distance_pair dt = this->sys->graph.ClosestAtom(*so2);
  printf ("the closest atom to the molecule is a distance of % 8.3f and is -- :\n", dt.first);
  dt.second->Print();

  for (Atom_it it = (*so2)->begin(); it != (*so2)->end(); it++) {
	(*it)->Print();
	dt = this->sys->graph.ClosestAtom(*it);
	printf ("closest atom (% 8.3f) -- \n", dt.first);
	dt.second->Print();
  }
  exit(1);

  /*
  for (Mol_it mol = sys_mols.begin(); mol != sys_mols.end(); mol++) {
	if (so2 == mol) continue;
	distances.push_back (std::make_pair (this->sys->Distance(*so2,*mol), *mol));
  }

  md_utility::pair_sort_first (distances.begin(), distances.end());
  //min_distances.push_back(distances[0]);

  fprintf (output, "% 8.3f % 8d %8s\n", distances[0].first, distances[0].second->MolID(), distances[0].second->Name().c_str());
  */

  return;

}

void Tester::DataOutput (const unsigned int timestep) { 
/*
  rewind(output);

  for (std::vector<mol_distance>::const_iterator it = min_distances.begin(); it != min_distances.end(); it++) {
	//std::cout << it->second->MolID() << std::endl;
	fprintf (output, "% 8.3f % 8d %s\n", it->first, it->second->MolID(), it->second->Name().c_str());
  }
  */

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

