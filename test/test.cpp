#include "test.h"


Tester::Tester (WaterSystemParams& wsp)
  :	
	Analyzer<XYZSystem> (wsp)
{
  return;
}

void Tester::Setup () {

  if (!timestep)
	rewind(output);

  this->sys->SetReparseLimit(1);

  s = Atom::FindByElement (sys_atoms, Atom::S);
  bondgraph::distance_vec Os = this->sys->graph.ClosestAtoms (s, 2, Atom::O, true);
  o1 = Os[0]->second;
  o2 = Os[1]->second;

  return;
}

void Tester::Analysis () {

  LoadAll();


  bondgraph::distance_vec closest_1 = this->sys->graph.ClosestAtoms (o1, 3);
  bondgraph::distance_vec closest_2 = this->sys->graph.ClosestAtoms (o2, 3);

  // data output 
  fprintf (output, "% 8d ", timestep);
  // atoms closest to o1
  fprintf (output, "o1 -- ");
  for (bondgraph::distance_vec::const_iterator it = closest_1.begin(); it != closest_1.end(); it++) {
	printf ("% 8.4f (%s/%d) ", it->first, it->second->Name().c_str(), it->second->ID());
  }

  // atoms closest to o2
  fprintf (output, "o2 -- ");
  for (bondgraph::distance_vec::const_iterator it = closest_2.begin(); it != closest_2.end(); it++) {
	printf ("% 8.4f (%s/%d) ", it->first, it->second->Name().c_str(), it->second->ID());
  }
  printf ("\n");

  return;
}

void Tester::DataOutput () {
  /*
	 rewind(output);

	 int i = 0;
	 for (bondgraph::distance_vec::const_iterator iti = so1.begin(), itj = so2.begin(); iti != so1.end(); iti++, itj++) {
	 fprintf (output, "%8d % 8.3f % 8.3f\n", ++i, iti->first, itj->first);
	 }

	 return; 
   */
}

void Tester::PostAnalysis () {

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

