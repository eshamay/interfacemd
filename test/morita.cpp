#include "morita.h"

SFGAnalyzer::SFGAnalyzer (int const argc, const char **argv, const WaterSystemParams& params)
	:	WaterSystem (argc, argv, params),
		sfg (SFGCalculator(&matrix)),
		MolChi (Complex_vec (0, complex<double>(0.0,0.0))),
		TimestepChi (Complex_vec (0, complex<double>(0.0,0.0))),
		TotalChi (Complex_vec (0, complex<double>(0.0,0.0)))
{

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

		// first let's find all the molecules in the interface
		this->FindWaters ();
		this->SliceWaters (45.0, 65.0);

		// and then update our bond data to reflect the interfacial region and find all the hydrogen bonds
		UpdateMatrix ();

		numMolsProcessed += int_mols.size();

		for (int mol = 0; mol < int_mols.size(); mol++) {

			water = static_cast<Water *>(int_mols[mol]);

			//int coord = static_cast<int>(matrix.WaterCoordination (water));
			//		printf ("%d\n", coord);

			MolChi.clear();

			// and then calculate the chi spectrum for the molecule SPS
			MolChi = sfg.Chi (*water);

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

		sys.LoadNext();
	}

	// final output of the data to the file
	OutputData ();

	fclose(output);

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
			double r = real(TotalChi[i]);
			double im = imag(TotalChi[i]);
			fprintf (output, "% 12.8e\t% 12.8e\t% 12.8e\n", freq, r/numMolsProcessed, im/numMolsProcessed);
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

	params.prmtop = "prmtop";
	params.mdcrd = "mdcrd";
	params.mdvel = "mdvel";
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
		params.output = "sfg.morita.dat";
	#endif
	params.posres = 0.100;
	params.pbcflip = 20.0;
	params.output_freq = 50;

	SFGAnalyzer analyzer (argc, argv, params);

	analyzer.Analyze ();

return 0;
}

