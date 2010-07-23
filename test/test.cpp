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

  //Mol_it so2 = Molecule::FindByType(sys_mols.begin(), sys_mols.end(), Molecule::SO2);
  Atom_it S = Atom::FindByElement (sys_atoms.begin(), sys_atoms.end(), Atom::S);
  bondgraph::distance_vec Os = this->sys->graph.ClosestAtoms (*S, 2, Atom::O, true);
  min_distances.push_back (Os[0]);
  min_distances.push_back (Os[1]);

  /*
  bondgraph::distance_vec dv = this->sys->graph.ClosestAtoms(*so2,3);
  for (bondgraph::distance_vec::const_iterator it = dv.begin(); it != dv.end(); it++) {
	printf ("% 8.3f % 8d %8s/%s\n", it->first, it->second->MolID(), it->second->Name().c_str(), it->second->Residue().c_str());
  }
  */

  //min_distances.push_back(dt);


  return;

}

void Tester::DataOutput () {
  rewind(output);

  for (bondgraph::distance_vec::const_iterator it = min_distances.begin(); it != min_distances.end(); it++) {
	fprintf (output, "% 8.3f % 8d %8s/%s\n", it->first, it->second->MolID(), it->second->Name().c_str(), it->second->Residue().c_str());
  }

  return; 
}

void Tester::PostAnalysis () {

  std::vector<double> dist;
  for (bondgraph::distance_vec::const_iterator it = min_distances.begin(); it != min_distances.end(); it++) {
	dist.push_back(it->first);
  }
  typedef std::vector< std::pair<double, int> > dist_t;
  dist_t histo = histogram::Histogram(dist.begin(), dist.end(), 200);

  rewind (output);

  for (dist_t::const_iterator it = histo.begin(); it != histo.end(); it++) {
	fprintf (output, "% 8.4f\t% 12d\n", it->first, it->second);
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

