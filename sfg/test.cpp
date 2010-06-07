#include "test.h"

Test::Test (WaterSystemParams& wsp)
  :	Analyzer<AmberSystem> (wsp)
{

  this->sys = new AmberSystem("prmtop", "mdcrd", "mdvel");

  return;
}

void Test::Setup () {

  this->LoadAll();
  this->LoadWaters();
  
  printf ("\nfirst\n");
  printf ("%d / %d atoms\n", (int)int_atoms.size(), (int)sys_atoms.size());
  printf ("%d / %d mols\n", (int)int_wats.size(), (int)sys_mols.size());

  std::pair<double,double> extents = this->ExtentPair();
  this->SliceWaters(int_wats, extents);
  this->UpdateAtoms (int_wats, int_atoms);

  printf ("second\n");
  printf ("%d / %d atoms\n", (int)int_atoms.size(), (int)sys_atoms.size());
  printf ("%d / %d mols\n", (int)int_wats.size(), (int)sys_mols.size());
  //printf ("atoms named: %s\n", int_atoms[0]->Name().c_str());

  return;
}

int main () {

  libconfig::Config cfg;
  cfg.readFile("system.cfg");

  std::string filename = cfg.lookup("analysis.test.filename");
  libconfig::Setting &analysis = cfg.lookup("analysis");
  analysis.add("filename", libconfig::Setting::TypeString) = filename;

  WaterSystemParams wsp (cfg);

  Test t (wsp);

  t.SystemAnalysis ();

  return 0;
}

