#include "densitytest.h"

// the command line arguments are the list of atom *NAMES* that are to be analyzed... exactly as they are named in the prmtop file
// the files prmtop, mdcrd, and mdvel should point to their respective files.
DensityAnalyzer::DensityAnalyzer (int const argc, const char **argv, const WaterSystemParams& params)
	:	WaterSystem (argc, argv, params)
{

	printf ("\n*** Performing a DENSITY analysis of the system ***\n");

	// from the command line, grab all the atom name that we'll be working with (remember that the first two arguments are interface locations
	atomNames.clear();
	#ifdef AVG
	for (int i = 3; i < argc; i++) {
	#else
	for (int i = 1; i < argc; i++) {
	#endif
		atomNames.push_back (argv[i]);
	}
	
	// some info before starting
	printf ("Analyzing densities of the atoms: ");
	RUN (atomNames) {
		printf ("%s  ", atomNames[i].c_str());
	}
	printf ("\n");

	// now we set up the histogram(s) such that we may bin the positions of each atom
	histo.resize(atomNames.size(), vector<int> (posbins, 0));
	
return;
}

void DensityAnalyzer::OutputData () {

	// starting from the beginning of the file (i.e. overwrite it)
	rewind (output);

	// the output value of the number density will be converted to the actual density in g/mL.
	// the volume of a slice of the slab
	double volume = Atom::Size()[x] * Atom::Size()[y] * Atom::Size()[z] / Atom::Size()[axis];
	volume *= posres;

	for (int i=0; i < posbins; i++) {
		fprintf (output, "% 10.4f", double(i)*posres+posmin);		// the bin's position value

		RUN2 (atomNames) {
			// The density is thus transformed. The resulting values are the densities (mol/mL) of each species.
			// The value of 1.6611293 comes from a combo of avogadro's number and the conversion to mL from angstroms^3.
			// To get the density in g/mL, just multiply this value by the molecular/atomic weight of the species.
			// For molarity (mol/L) multiply this value by 1000 (to convert from mL to L)
			double density = double(histo[j][i]) / volume / ((double)timestep+1.0) * 1.6611293;
			#ifdef AVG
			density /= 2.0;		// if we're averaging 2 interfaces, then these values need to be halved
			#endif

			fprintf (output, "% 13.7f", density);			// for each atom printout the histogram value at that position
		}

		fprintf (output, "\n");
	}

return;
}

vector<int> DensityAnalyzer::AtomDensity (string const atomname) {

	vector<int> density (posbins, 0);
	
	// and now run through the actual calculations to find number densities
	RUN (sys) {
			
		// find the atom of interest
		Atom * patom = sys[i];
		
		//if (_sys[i]->Name().find(atomname) == string::npos) continue;
		if (sys[i]->Name() != atomname) continue;

		// grab the position info of the atom
		VecR r = patom->Position();
		double position = r[axis];
		if (position < pbcflip) position += Atom::Size()[axis];

		#ifdef AVG
		// here the bin will be selected based on the distance to a given interface. Negative distances are inside the water phase, positive are in the CCl4
		double distance = (position > middle) ? position - int_high : int_low - position;
		int bin = (int)((distance - posmin)/posres);
		#else
		//if (position < START or position > END) continue;		// only bin stuff within the bounds that have been set up
		int bin = (int)((position - posmin)/posres);
		#endif

		if (bin > posbins) {
			printf ("DensityAnalyzer::AtomDensity - something is wrong - not enough bins in the histogram for analysis?\n");
			exit(1);
		}
		// and bin it into the density histogram
		density[bin]++;

	}

return (density);
}

void DensityAnalyzer::SystemDensities () {

	// now let's run through the timestesp
	for (timestep=0; timestep < timesteps; timestep++) {
		
		// for each atom that we're testing we'll add the histogram data into the final data-set
		for (unsigned int atom = 0; atom < atomNames.size(); atom++) {
			vector<int> atomDensity = this->AtomDensity (atomNames[atom]);

			// once we have the data for each atom for each timestep, let's add it into the running total
			for (unsigned int bin = 0; bin < posbins; ++bin) {
				histo[atom][bin] += atomDensity[bin];
			}
		}
			
		// and set up the system for the next timestep
		sys.LoadNext();

		// output some info
		this->OutputStatus ();
		this->OutputData ();

	}

	this->OutputData ();

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
	params.mdvel = "";
	params.axis = y;
	params.timesteps = 200000;
	#ifdef RESTART
		params.restart = 138750;
	#else
		params.restart = 0;
	#endif

	#ifdef AVG
		params.avg = true;
		params.posmin = -40.0;
		params.posmax = 40.0;
		params.output = "density.avg.dat";
	#else
		params.avg = false;
		params.posmin = -5.0;
		params.posmax = 150.0;
		params.output = "density.dat";
	#endif
	params.posres = 0.100;
	params.pbcflip = 15.0;
	params.output_freq = 25;


	DensityAnalyzer den (argc, argv, params);

	den.SystemDensities();

return 0;
}

