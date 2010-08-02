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

  LoadAll();
  this->sys->SetReparseLimit(1);

  s = Atom::FindByElement (sys_atoms, Atom::S);
  bondgraph::distance_vec Os = this->sys->graph.ClosestAtoms (s, 2, Atom::O, true);
  o1 = Os[0].second;
  o2 = Os[1].second;

  return;
}

void Tester::Analysis () {
  LoadAll();

  /*
  // find out the hbond-distance factor for each O atom. 
  // The Distance Factor is a dimensionless number that represents the position of a hydrogen between two oxygens.
  // An H exactly in the middle of 2 oxygens will have Q=0.0
  // an H sitting right on top of the water (target) oxygen will have a value of 1.0, and sitting right on top of the source (SO2) oxygen will have a value of -1.0
  // Q = (dOH1 - dOH2)/D
  // dOH1 = distance of H to SO2 Oxygen
  // dOH2 = distance of H to H2O oxygen
  // D = total distance
  bondgraph::distance_pair h1_pair = this->sys->graph.ClosestAtom (o1, Atom::H);
  bondgraph::distance_pair h2_pair = this->sys->graph.ClosestAtom (o2, Atom::H);

  AtomPtr target_o1 = h1_pair.second->ParentMolecule()->GetAtom(Atom::O);
  AtomPtr target_o2 = h2_pair.second->ParentMolecule()->GetAtom(Atom::O);

  double target_distance1 = this->sys->graph.Distance(target_o1, h1_pair.second);
  double target_distance2 = this->sys->graph.Distance(target_o2, h2_pair.second);

  double Q1 = (h1_pair.first - target_distance1)/(h1_pair.first + target_distance1);
  double Q2 = (h2_pair.first - target_distance2)/(h2_pair.first + target_distance2);

  // data output 
  fprintf (output, "% 8d % 8.4f % 8.4f\n", timestep, Q1, Q2);
  */

  // atoms closest to o1
  fprintf (output, "o1 -- ");
  for (bondgraph::distance_vec::const_iterator it = closest_1.begin(); it != closest_1.end(); it++) {
	fprintf (output, "% 8.4f (%s/%d) ", it->first, it->second->Name().c_str(), it->second->ID());
  }

  // atoms closest to o2
  fprintf (output, "o2 -- ");
  for (bondgraph::distance_vec::const_iterator it = closest_2.begin(); it != closest_2.end(); it++) {
	fprintf (output, "% 8.4f (%s/%d) ", it->first, it->second->Name().c_str(), it->second->ID());
  }
  fprintf (output, "\n");

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

