#include "morita.h"

SFGCalculator SFGAnalyzer::sfg;
Complex_vec SFGAnalyzer::TimestepChi;
Complex_vec SFGAnalyzer::TotalChi;

unsigned long SFGAnalyzer::numMolsProcessed = 0;
bool SFGAnalyzer::firstMol = true;
bool SFGAnalyzer::firstTimeStep = true;


SFGAnalyzer::SFGAnalyzer (WaterSystemParams& wsp)
  :	Analyzer<AmberSystem> (wsp)
	//sfg (SFGCalculator(&this->_graph))
{

  //this->sys = new AmberSystem("prmtop", "mdcrd", "mdvel");

  printf ("\n*** Performing an SFG analysis of the system ***\n");

  return;
}

void SFGAnalyzer::Setup () {

  // Load up all the water molecules and atoms
  this->LoadWaters();
  int numWats = (int)int_wats.size();
  //printf ("%d/%d molecules are waters\n", numWats, (int)sys_mols.size());

  // keep only the waters within a given region of the slab for analysis
  std::pair<double,double> extents = std::make_pair<double,double> (
	  WaterSystem<AmberSystem>::posmin,
	  WaterSystem<AmberSystem>::posmax
	  );
  this->SliceWaters(int_wats, extents);
  printf ("%d/%d waters, %d/%d atoms\n", (int)int_wats.size(), numWats, (int)int_atoms.size(), (int)sys_atoms.size());

  return;
}

void SFGAnalyzer::Analysis () {

  TimestepChi.clear();	// it's a new timestep
  firstMol = true;		// every timestep we will have to go through all the molecules again

  // and then update our bond data to reflect the interfacial region and find all the hydrogen bonds
  UpdateGraph ();

  std::for_each(int_wats.begin(), int_wats.end(), this->SFGProcess);

  // we collect the data for each timestep into the running total for the end spectrum
  CollectChi (TotalChi.begin(), TotalChi.end(), TimestepChi.begin());


  return;
}

// take each of the interfacial waters and flip the molecules such that they are mirrored about a plane that runs through the oxygen and is normal to the given axis.
/*
   void SFGAnalyzer::FlipWaters (const coord axis) {
   RUN (int_wats) {
   int_wats[i]->Flip(axis);
   }
   return;
   }
   */

// output data to the file
void SFGAnalyzer::DataOutput (const unsigned int timestep) {

  rewind (output);

  double freq = START_FREQ * AU2WAVENUMBER;
  double scale = (double)numMolsProcessed;
  double r, im;

  for (Complex_vec::const_iterator it = TotalChi.begin(); it != TotalChi.end(); it++) {
	r = real(*it) / scale;
	im = imag(*it) / scale;

	fprintf (output, "% 20.4e% 20.4e% 20.4e\n", freq, r, im);

	freq += FREQ_STEP * AU2WAVENUMBER;
  }

  fflush (output);

  // reload the interface waters once in a while because they tend to move in and out of the interfacial region
  this->Setup();

  return;
}

void SFGAnalyzer::PostAnalysis ()
{
  return;
}

int main () {

  libconfig::Config cfg;
  cfg.readFile("system.cfg");

  std::string filename = cfg.lookup("analysis.sfg.filename");
  libconfig::Setting &analysis = cfg.lookup("analysis");
  analysis.add("filename", libconfig::Setting::TypeString) = filename;

  WaterSystemParams wsp (cfg);

  SFGAnalyzer analyzer (wsp);

  analyzer.SystemAnalysis ();

  return 0;
}

