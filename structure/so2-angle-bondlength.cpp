#include "so2-angle-bondlength.h"


Tester::Tester (WaterSystemParams& wsp)
  :	
	Analyzer<XYZSystem> (wsp)
{
  return;
}

void Tester::Setup () {

  LoadAll();
  this->sys->SetReparseLimit(1);

  return;
}

void Tester::Analysis () {
  LoadAll();

  MolPtr mol = Molecule::FindByType(sys_mols, Molecule::SO2);
  so2 = new SulfurDioxide(mol);
  so2->SetAtoms();

  so1_data.push_back(so2->SO1().norm());
  so2_data.push_back(so2->SO2().norm());
  double theta = so2->SO1() < so2->SO2();
  theta_data.push_back(acos(theta)*180.0/M_PI);
  delete so2;

 return;
}

void Tester::DataOutput () {
  rewind(output);

  fprintf (output, "timestep so1length so2length angle\n");
  for (int i = 0; i < so1_data.size(); i++) {
	fprintf (output, "% 8d % 8.4f % 8.4f % 8.4f\n", i+1, so1_data[i], so2_data[i], theta_data[i]);
  }

}

void Tester::PostAnalysis () {
  return;
}


int main () {
  libconfig::Config cfg;
  cfg.readFile("system.cfg");

  std::string filename ("so2-bond-angle.dat");
  libconfig::Setting &analysis = cfg.lookup("analysis");
  analysis.add("filename", libconfig::Setting::TypeString) = filename;

  WaterSystemParams wsp (cfg);

  Tester test (wsp);

  test.SystemAnalysis ();

  return 0;
}

