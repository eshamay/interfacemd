#include "morita.h"

SFGAnalyzer::SFGAnalyzer (int const argc, const char **argv, const WaterSystemParams& params)
	:	WaterSystem<AmberSystem> (params),
		sfg (SFGCalculator(&this->graph)),
		MolChi (Complex_vec (0, complex<double>(0.0,0.0))),
		TimestepChi (Complex_vec (0, complex<double>(0.0,0.0))),
		TotalChi (Complex_vec (0, complex<double>(0.0,0.0)))
{

	this->sys = new AmberSystem("prmtop", "mdcrd", "mdvel");

	printf ("\n*** Performing an SFG analysis of the system ***\n");

return;
}

void SFGAnalyzer::Analyze () {

	Water * water;	// our prototypical water molecule

	bool firstmol = true;			// first molecule processed

	numMolsProcessed = 0;

	// okay, now let's start iterating over timesteps
	for (timestep = 0; timestep < timesteps; timestep++) {

		TimestepChi.clear();	// it's a new timestep
		firstmol = true;		// every timestep we will have to go through all the molecules again

		RUN (sys->Molecules()) {
			//Molecule * mol = sys.Molecules(i);
			//mol->Print();
		}

		// first let's find all the molecules in the interface
		this->FindWaters ();
		this->SliceWaters (55.0, 75.0);

		// and then update our bond data to reflect the interfacial region and find all the hydrogen bonds
		UpdateGraph ();
		// only grab the OH-waters for now
		//this->SliceWaterCoordination (OOH);

		for (int mol = 0; mol < (int)int_mols.size(); mol++) {

			water = static_cast<Water *>(int_mols[mol]);

			MolChi.clear();

			// and then calculate the chi spectrum for the molecule SPS
			MolChi = sfg.Beta (*water, 0,2,1);

			numMolsProcessed++;

			// when starting a new timestep...
			if (firstmol) {
				TimestepChi.resize (MolChi.size(), complex<double> (0.0, 0.0));
				firstmol = false;
			}

			// for the very first timestep...
			if (!timestep) {
				int numData = MolChi.size();	// number of data points collected for the chi spectrum
				TotalChi.resize (numData, complex<double> (0.0, 0.0));
			}

			// perform the summation for averaging over the system
			CollectChi (MolChi, TimestepChi);
		}

		// we collect the data for each timestep into the running total
		CollectChi (TimestepChi, TotalChi);

		// now output something to the screen (once in a while = every set # of timesteps)
		OutputStatus ();

		// and once in a while, also output the data to a file for reading
		OutputData ();

		sys->LoadNext();
	}

	// final output of the data to the file
	OutputData ();

return;
}

// sum up the Chi spectra from each molecule
void SFGAnalyzer::CollectChi (Complex_vec& newchi, Complex_vec& totalchi) {

	RUN (newchi) {
		totalchi[i] += newchi[i];
	}

return;
}

// output data to the file
void SFGAnalyzer::OutputData () {

	if (timestep && !(timestep % (output_freq * 10))) {
		rewind (output);

		//RUN (TotalChi) {
			//printf ("t - %f\n", abs(TotalChi[i]));
		//}
		RUN (TotalChi) {
			double freq = (double(i)*FREQ_STEP+START_FREQ)*AU2WAVENUMBER;
			double scale = (double)numMolsProcessed;
			double r = real(TotalChi[i]) / scale;
			double im = imag(TotalChi[i]) / scale;
			fprintf (output, "% 20.8e\t% 20.8e\t% 20.8e\n", freq, r, im);
		}

		fflush (output);
	}

return;
}

int main (const int argc, const char **argv) {

	#ifdef AVG
		if (argc < 3) {
			printf ("for averaging, we need the two interface locations: <low> <high>\n");
			exit(1);
		}
	#endif

	WaterSystemParams params;

	params.axis = y;
	params.timesteps = 200000;
	#ifdef RESTART
		params.restart = 100000;
	#endif
	#ifdef AVG
		params.avg = true;
		params.posmin = -40.0;
		params.posmax = 40.0;
		params.output = "sfg.morita.avg.dat";
	#else
		params.avg = false;
		params.posmin = -5.0;
		params.posmax = 150.0;
		params.output = "sfg.alpha-switch.dat";
	#endif
	params.posres = 0.100;
	params.pbcflip = 15.0;
	params.output_freq = 50;

	SFGAnalyzer analyzer (argc, argv, params);

	analyzer.Analyze ();

return 0;
}

