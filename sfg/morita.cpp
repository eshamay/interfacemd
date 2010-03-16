#include "morita.h"

SFGAnalyzer::SFGAnalyzer (WaterSystemParams& params)
  :	Analyzer (params),
	sfg (SFGCalculator(&this->_graph)),
	Molecular_Beta (Complex_vec (0, complex<double>(0.0,0.0))),
	TimestepChi (Complex_vec (0, complex<double>(0.0,0.0))),
	TotalChi (Complex_vec (0, complex<double>(0.0,0.0))),
	numMolsProcessed(0), firstmol(true), firsttimestep(true)
{

  this->sys = new AmberSystem("prmtop", "mdcrd", "mdvel");

  printf ("\n*** Performing an SFG analysis of the system ***\n");

  return;
}

void SFGAnalyzer::Setup () {
  return;
}

void SFGAnalyzer::Analysis () {


  TimestepChi.clear();	// it's a new timestep
  firstmol = true;		// every timestep we will have to go through all the molecules again

  // first find all the waters in the system
  FindWaters();
  // first let's find all the waters in the interface
  SLICE_BY_POSITION(int_wats, Water *, 60.0, 80.0);


  // and then update our bond data to reflect the interfacial region and find all the hydrogen bonds
  UpdateGraph ();

  // only grab the OH-waters for now
  //this->SliceWaterCoordination (OOH);

  Water * water;
  for (int mol = 0; mol < (int)int_wats.size(); mol++) {

	water = int_wats[mol];

	Molecular_Beta.clear();

	// and then calculate the chi spectrum for the molecule 
	// 0,2,1 = SSP. S = X and Z axes, P = Y axis
	Molecular_Beta = sfg.Beta (*water, 0,2,1);

	numMolsProcessed++;

	// when starting a new timestep...
	if (firstmol) {
	  TimestepChi.resize (Molecular_Beta.size(), complex<double> (0.0, 0.0));
	  firstmol = false;
	}

	// for the very first timestep...
	if (firsttimestep) {
	  TotalChi.resize (Molecular_Beta.size(), complex<double> (0.0, 0.0));
	  firsttimestep = false;
	}

	// perform the summation for averaging over the system
	CollectChi (Molecular_Beta, TimestepChi);
  }

  // we collect the data for each timestep into the running total
  CollectChi (TimestepChi, TotalChi);

  return;
}

// sum up the Chi spectra from each molecule
void SFGAnalyzer::CollectChi (Complex_vec& newchi, Complex_vec& totalchi) {

  RUN (newchi) {
	totalchi[i] += newchi[i];
  }

  return;
}


// take each of the interfacial waters and flip the molecules such that they are mirrored about a plane that runs through the oxygen and is normal to the given axis.
void SFGAnalyzer::FlipWaters (const coord axis) {
  RUN (int_wats) {
	int_wats[i]->Flip(axis);
  }
  return;
}

// output data to the file
void SFGAnalyzer::DataOutput (const unsigned int timestep) {

  if (timestep && !(timestep % (output_freq * 10))) {
	rewind (output);

	RUN (TotalChi) {
	  double freq = (double(i)*FREQ_STEP+START_FREQ)*AU2WAVENUMBER;
	  double scale = (double)numMolsProcessed;
	  double r = real(TotalChi[i]) / scale;
	  double im = imag(TotalChi[i]) / scale;
	  fprintf (output, "% 20.4e% 20.4e% 20.4e\n", freq, r, im);
	}

	fflush (output);
  }

  return;
}

void SFGAnalyzer::PostAnalysis ()
{
  return;
}

int main () {


  WaterSystemParams params ("Morita-SFG.dat", 100000);

  SFGAnalyzer analyzer (params);

  analyzer.SystemAnalysis ();

  return 0;
}

