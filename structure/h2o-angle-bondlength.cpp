#include "h2o-angle-bondlength.h"


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
  LoadWaters();

  for (Mol_it it = int_wats.begin(); it != int_wats.end(); it++) {
	wat = new Water(*it);
	wat->SetAtoms();

	oh.push_back(wat->OH1()->norm());
	oh.push_back(wat->OH2()->norm());
	double angle = wat->Angle();
	theta.push_back(acos(angle)*180.0/M_PI);
	delete wat;
  }

  return;
}

void Tester::DataOutput () {
  rewind(output);

  histogram::histogram_t oh_histo = histogram::Histogram (oh.begin(), oh.end(), 200);
  int oh_histo_max = histogram::MaxPopulation (oh_histo.begin(), oh_histo.end());
  histogram::histogram_t theta_histo = histogram::Histogram (theta.begin(), theta.end(), 200);
  int theta_histo_max = histogram::MaxPopulation (theta_histo.begin(), theta_histo.end());

  for (unsigned i = 0; i < oh_histo.size(); i++) {
	fprintf (output, "% 8.4f % 8.4f % 8.4f % 8.4f\n", 
		oh_histo[i].first, (double)oh_histo[i].second/oh_histo_max, 
		theta_histo[i].first, (double)theta_histo[i].second/theta_histo_max);
  }

}

void Tester::PostAnalysis () {
  return;
}


int main () {
  libconfig::Config cfg;
  cfg.readFile("system.cfg");

  std::string filename ("h2o-bond-angle.dat");
  libconfig::Setting &analysis = cfg.lookup("analysis");
  analysis.add("filename", libconfig::Setting::TypeString) = filename;

  WaterSystemParams wsp (cfg);

  Tester test (wsp);

  test.SystemAnalysis ();

  return 0;
}

