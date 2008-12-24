#include "morita.h"

int main (int argc, char **argv) {

	// here is our system for analysis
	AmberSystem sys (PRMTOP, MDCRD, FORCE);

	// This will be our sfg-analysis toybox
	SFGWaterAnalyzer sfg;

	// a little info sent to the screen about what we're doing
	OutputHeader ();

	FILE * output = (FILE *)NULL;

	output = fopen (OUTPUTFILE, "w");
	if (output == (FILE *)NULL) {
		std::cout << "couldn't open the output file!!" << endl;
		exit (1);
	}

	// let's set up a few useful containers that we'll use later on
	std::vector< complex<double> > MolChi;			// chi for a given molecule
	std::vector< complex<double> > TimestepChi;		// chi spectrum of several different molecules collected over an entire timestep
	std::vector< complex<double> > TotalChi;		// Running total of all the data for several timesteps

	Water * water;	// our prototypical water molecule

	// Sets up the axes for the system
	const int P = 2;		// perpendicular to the interface
	const int S1 = 0;		// the two that are parallel to the surface
	const int S2 = 1;


	// just so we don't have to always iterate over the entire system, we're only going to look at interfacial molecules - the ones were interested in.
	std::vector<Water *> int_mols;
	// and these are the atoms of those molecules
	std::vector<Atom *> int_atoms;
	
	bool firsttimestep = true;	// first time through (first timestep)
	bool firstmol;			// first molecule processed

	// okay, now let's start iterating over timesteps
	for (int step=0; step < TIMESTEPS; step++) {

		TimestepChi.clear();	// it's a new timestep
		firstmol = true;		// every timestep we will have to go through all the molecules again

		// first let's find all the molecules in the interface
		FindInterfacialAtoms (int_mols, int_atoms, sys);

		// and then update our bond data to reflect the interfacial region and find all the hydrogen bonds
		sys.UpdateGraph (int_atoms);

		for (int mol = 0; mol < int_mols.size(); mol++) {
			
			water = int_mols[mol];

			sfg.Reset();		// reset the calculator for a new molecule
			MolChi.clear();

			// and then calculate the chi spectrum for the molecule SPS
			MolChi = sfg.Chi (*water, S1, P, S1);

			// when starting a new timestep...
			if (firstmol) {
				TimestepChi.resize (MolChi.size(), complex<double> (0.0, 0.0));
				firstmol = false;
			}

			// for the very first timestep...
			if (firsttimestep) {
				int numData = MolChi.size();	// number of data points collected for the chi spectrum
				TotalChi.resize (numData, complex<double> (0.0, 0.0));
				firsttimestep = false;
			}
		
			// perform the summation for averaging over the system
			CollectChi (MolChi, TimestepChi);

			// now add in the next equivalent polarization
			MolChi = sfg.Chi (*water, S2, P, S2);
	
			// perform the summation for the other half of the polarization
			CollectChi (MolChi, TimestepChi);
		}
		
		// now output something to the screen (once in a while = every set # of timesteps)
		OutputStatus (step);

		// we collect the data for each timestep into the running total
		CollectChi (TimestepChi, TotalChi);

		// and once in a while, also output the data to a file for reading
		if (!(step % (OUTPUT_FREQ * 10))) {
			OutputData (output, TotalChi);		// output the data to the data file
		}

		sys.LoadNext();
	}

	// final output of the data to the file
	OutputData (output, TotalChi);

	fclose(output);

return 0;
}

// sum up the Chi spectra from each molecule
void CollectChi (std::vector< std::complex<double> >& Chi_step, std::vector< std::complex<double> >& Chi_total) {

	RUN (Chi_step) {
		Chi_total[i] += Chi_step[i];
	}

return;
}
	
// print an informative header
void OutputHeader () {
	
	printf ("Generating SFG for %d timesteps\n\"*\" = %d steps\n", TIMESTEPS, OUTPUT_FREQ);
	printf ("Polarization = %s\n", POLARIZATION);
	printf ("outputting to file: %s\n", OUTPUTFILE);

return;
}

// output a status meter to the screen
void OutputStatus (int const count) {

	if (!(count % (OUTPUT_FREQ * 10)))
		printf ("\n%10d/%d)  ", count, TIMESTEPS);

	if (!(count % OUTPUT_FREQ)) 
		printf ("*");
	
	fflush (stdout);


return;
}

// output data to the file
void OutputData (FILE * fp, vector< complex<double> >& chi) {

	rewind (fp);

	RUN (chi) {
		double freq = (double(i)*FREQ_STEP+START_FREQ)*FREQFACTOR;
		fprintf (fp, "% 12.8e\t% 12.8e\n", freq, abs(chi[i])*abs(chi[i]));
	}

	fflush (fp);

return;
}

void FindInterfacialAtoms (vector<Water *>& int_mols, vector<Atom *>& int_atoms, AmberSystem& sys) {
	
	int_mols.clear();
	int_atoms.clear();

	Molecule * pmol;

	// go through the system
	RUN (sys.Molecules()) {
		// grab each molecule
		pmol = sys.Molecules(i);

		// we're only looking at waters for SFG analysis right now
		if (pmol->Name() != "h2o") continue;

		// first thing we do is grab a water molecule to work with
		Water * water = static_cast<Water *>(pmol);

		double position = water->Atoms(0)->Position()[axis];
		// and find molecules that sit within the interface.
		if (position < PBC_FLIP) position += Atom::Size()[axis];		// adjust for funky boundaries
		// these values have to be adjusted for each system
		if (position < INTERFACE_LOW or position > INTERFACE_HIGH) continue;				// set the interface cutoffs
		
		int_mols.push_back (water);
		RUN2(water->Atoms()) {
			int_atoms.push_back (water->Atoms(j));
		}
	}

return;
}
