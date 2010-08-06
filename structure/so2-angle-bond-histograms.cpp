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

  so.push_back(so2->SO1().norm());
  so.push_back(so2->SO2().norm());
  double angle = so2->Angle();
  theta.push_back(acos(angle)*180.0/M_PI);
  delete so2;

  return;
}

void Tester::DataOutput () {
  rewind(output);

  histogram::histogram_t so_histo = histogram::Histogram (so.begin(), so.end(), 100);
  int so_histo_max = histogram::MaxPopulation (so_histo.begin(), so_histo.end());
  histogram::histogram_t theta_histo = histogram::Histogram (theta.begin(), theta.end(), 100);
  int theta_histo_max = histogram::MaxPopulation (theta_histo.begin(), theta_histo.end());

  for (unsigned i = 0; i < so_histo.size(); i++) {
	fprintf (output, "% 8.4f % 8.4f % 8.4f % 8.4f\n", 
		so_histo[i].first, (double)so_histo[i].second/so_histo_max, 
		theta_histo[i].first, (double)theta_histo[i].second/theta_histo_max);
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

