#define MPI_SYS
#include "morita.h"

int main (int argc, char **argv) {

	// here is our system for analysis
	MPIMolSystem sys (&argc, &argv, PRMTOP, MDCRD, FORCE);

	// This will be our sfg-analysis toybox
	SFGWaterAnalyzer sfg;

	// a little info sent to the screen about what we're doing
	OutputHeader (sys);

	FILE * output = (FILE *)NULL;

  if (sys.Master())	{	// only one mpi node needs to open a file for writing data
	output = fopen (OUTPUTFILE, "w");
	if (output == (FILE *)NULL) {
		std::cout << "couldn't open the output file!!" << endl;
		exit (1);
	}
  }

	// let's set up a few useful containers that we'll use later on
	std::vector< complex<double> > MolChi;			// chi for a given molecule
	std::vector< complex<double> > TimestepChi;		// chi spectrum of several different molecules collected over an entire timestep
	std::vector< complex<double> > TotalChi;		// Running total of all the data for several timesteps

	Water * water;	// our prototypical water molecule

	// Sets up the axes for the system
	const int P = 1;		// perpendicular to the interface
	const int S1 = 0;		// the two that are parallel to the surface
	const int S2 = 2;


	// just so we don't have to always iterate over the entire system, we're only going to look at interfacial molecules - the ones were interested in.
	std::vector<Molecule *> int_wats;
	// and these are the atoms of those molecules
	std::vector<Atom *> int_atoms;

	bool firsttimestep = true;	// first time through (first timestep)
	bool firstmol;			// first molecule processed

	// okay, now let's start iterating over timesteps
	for (int step=0; step < TIMESTEPS; step++) {

		TimestepChi.clear();	// it's a new timestep
		firstmol = true;		// every timestep we will have to go through all the molecules again

		// first let's find all the molecules in the interface
		FindInterfacialAtoms (int_wats, int_atoms, sys);

/*
for (int i = 0; i < 5; i++) {
	printf ("%d) %d %d", sys.ID(), i, int_atoms[i]->ID());
}
*/

	  if (sys.Master())
		// and then update our bond data to reflect the interfacial region and find all the hydrogen bonds
		sys.UpdateGraph (int_atoms);

		// each node needs to know which molecules to analyze
		int start_mol = BLOCK_LOW(sys.ID(), sys.Procs(), int_wats.size());
		int end_mol = BLOCK_HIGH(sys.ID(), sys.Procs(), int_wats.size());
		//printf ("node %d of %d covers [%d - %d]\n", sys.ID(), sys.Procs(), start_mol, end_mol);

		for (int mol = start_mol; mol <= end_mol; mol++) {
			// we're only looking at waters for SFG analysis right now
			if (int_wats[mol]->Name() != "h2o") continue;

			// first thing we do is grab a water molecule to work with
			water = static_cast<Water *>(int_wats[mol]);

			sfg.Reset();		// reset the calculator for a new molecule
			MolChi.clear();

			// and then calculate the chi spectrum for the molecule
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
		OutputStatus (step, sys);

		// if doing MPI runs, then collect all the timestep data from each node
		MPI_PackageChi(TimestepChi, TotalChi, sys);

		// and once in a while, also output the data to a file for reading
		if (!(step % (OUTPUT_FREQ * 10))) {
		  if (sys.Master())
			OutputData (output, TotalChi);		// output the data to the data file
		}

		sys.LoadNext();
	}

	// final output of the data to the file
  if (sys.Master()) {
	OutputData (output, TotalChi);

	fclose(output);
  }

return 0;
}

// package up Chi from each mpi node such that it can be shipped off to the master for publication
void MPI_PackageChi (std::vector< std::complex<double> >& TimestepChi, std::vector< std::complex<double> >& TotalChi, MPIMolSystem& sys) {

	int numData = TimestepChi.size();

	double real_data[numData];
	double imag_data[numData];

	// when running mpi, each node needs to deconstruct the TimestepChi in order to ship the data to the master node after each timestep
	//  here we package each node's data for the timestep
	RUN (TimestepChi) {
		real_data[i] = real(TimestepChi[i]);
		imag_data[i] = imag(TimestepChi[i]);
	}

	double total_real [numData];
	double total_imag [numData];
	for (int i=0; i < numData; i++) {
		total_real[i] = 0.0;
		total_imag[i] = 0.0;
	}

	// now we sum the data from each node
	MPI_Allreduce (real_data, total_real, numData, MPI_DOUBLE, MPI_SUM, sys.WorldComm());
	MPI_Allreduce (imag_data, total_imag, numData, MPI_DOUBLE, MPI_SUM, sys.WorldComm());

	// and then the master takes the total sum for the timestep and adds it into the running total
	if (sys.Master()) {
		RUN(TotalChi) {
			TotalChi[i] += complex<double> (total_real[i], total_imag[i]);
		}
	}

return;
}

// sum up the Chi spectra from each molecule
void CollectChi (std::vector< std::complex<double> >& Chi_step, std::vector< std::complex<double> >& Chi_total) {

	RUN (Chi_step) {
		Chi_total[i] += Chi_step[i];
	}

return;
}

// print an informative header
void OutputHeader (MPIMolSystem& sys) {

  if (sys.Master()) {
	printf ("Generating SFG for %d timesteps\n\"*\" = %d steps\n", TIMESTEPS, OUTPUT_FREQ);
	printf ("Polarization = %s\n", POLARIZATION);
	printf ("outputting to file: %s\n", OUTPUTFILE);
  }

return;
}

// output a status meter to the screen
void OutputStatus (int const count, MPIMolSystem& sys) {

  if (sys.Master()) {
	if (!(count % (OUTPUT_FREQ * 10)))
		printf ("\n%10d/%d)  ", count, TIMESTEPS);

	if (!(count % OUTPUT_FREQ))
		printf ("*");

	fflush (stdout);

  }

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

void FindInterfacialAtoms (vector<Molecule *>& int_wats, vector<Atom *>& int_atoms, MPIMolSystem& sys) {

	int_wats.clear();
	int_atoms.clear();

	Molecule * pmol;

  if (sys.Master()) {
	// go through the system
	RUN (sys.Molecules()) {
		pmol = sys.Molecules(i);
		double position = pmol->Atoms(0)->Position()[axis];
		// and find molecules that sit within the interface.
		if (position < PBC_FLIP) position += Atom::Size()[axis];		// adjust for funky boundaries
		// these values have to be adjusted for each system
		if (position < INTERFACE_LOW || position > INTERFACE_HIGH) continue;				// set the interface cutoffs

		int_wats.push_back (pmol);
		RUN2(pmol->Atoms()) {
			int_atoms.push_back (pmol->Atoms(j));
		}
	}

  }
// now since each node has its own memory, we can't send Mol * or Atom *, so we send ID #s

	// first gather all the atom and molecule IDs
	vector<int> molIDs, atomIDs;
	if (sys.Master()) {

		RUN (int_atoms)
			atomIDs.push_back(int_atoms[i]->ID());
		RUN (int_wats)
			molIDs.push_back(int_wats[i]->MolID());
	}

	// deconstruct, ship out, and reconstruct on each node
	sys.BroadcastVector (atomIDs);
	sys.BroadcastVector (molIDs);

	// then the other nodes fix up their own lists of atoms and molecules in the interface
	if (!sys.Master()) {

		RUN (molIDs)
			int_wats.push_back (sys.Molecules (molIDs[i]));
		RUN (atomIDs)
			int_atoms.push_back (sys.Atoms (atomIDs[i]));
	}

return;
}
