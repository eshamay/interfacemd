#include "charge-histogram.h"

template <class U>
void ChargeBinner<U>::operator() (Atom * t) {
  _histogram[Analyzer<U>::PositionBin(t)] += AtomCharge(t);
  return;
}

template <class U>
double ChargeBinner<U>::AtomCharge (Atom * t) {
  double charge;

  if (t->Name().find("H") != std::string::npos) charge = 0.365;
  //else if (t->Name().find("O1") != std::string::npos) charge = -0.65;
  //else if (t->Name().find("O2") != std::string::npos) charge = -0.65;
  //else if (t->Name().find("O3") != std::string::npos) charge = -0.65;
  else if (t->Name().find("S") != std::string::npos) charge = 2.0;
  else if (t->Name().find("O1") != std::string::npos) charge = -1.0;
  else if (t->Name().find("O2") != std::string::npos) charge = -1.0;
  else if (t->Name().find("O3") != std::string::npos) charge = -1.0;
  else if (t->Name().find("O4") != std::string::npos) charge = -1.0;
  else if (t->Name().find("O") != std::string::npos) charge = -0.73;
  else if (t->Name().find("Cl1") != std::string::npos) charge = 0.0404;
  else if (t->Name().find("Cl2") != std::string::npos) charge = 0.0404;
  else if (t->Name().find("Cl3") != std::string::npos) charge = 0.0404;
  else if (t->Name().find("Cl4") != std::string::npos) charge = 0.0404;
  else if (t->Name().find("Cl") != std::string::npos) charge = -1.0;
  else if (t->Name().find("C") != std::string::npos) charge = -0.1616;
  else if (t->Name().find("NA") != std::string::npos) charge = 1.0;
  else if (t->Name().find("N") != std::string::npos) charge = 0.95;

  return charge;
}

template <class U>
// Output data from the histograms
void ChargeBinner<U>::Output (FILE * output, const int timestep) {

  rewind (output);

  double scale = (!timestep) ? 1.0 : double(timestep);

  // first print out a header consisting of column headers
  fprintf (output, "Position\tCharge\n");	// The position column

  double position;
  double charge;
  // for every position in the system
  for (int pos = 0; pos < Analyzer<U>::posbins; pos++) {
    // print the slab position
    position = double(pos)*Analyzer<U>::posres + Analyzer<U>::posmin;
    charge = _histogram[pos] / scale;
    fprintf (output, "% 13.5f% 13.5f\n", position, charge);
  }
return;
}

int main () {

  libconfig::Config cfg;
  cfg.readFile("system.cfg");

  std::string filename = cfg.lookup("analysis.charge.filename");
  libconfig::Setting &analysis = cfg.lookup("analysis");
  analysis.add("filename", libconfig::Setting::TypeString) = filename;

  WaterSystemParams wsp (cfg);

  ChargeAnalyzer<AmberSystem> ca (wsp);

  ca.SystemAnalysis();

  return 0;
}
