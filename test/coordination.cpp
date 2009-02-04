#include "coordination.h"

CoordinationTest::CoordinationTest () {
	
	// here is our system for analysis
	sys = new AmberSystem (PRMTOP, MDCRD, FORCE);

	printf ("Running a test to find water coordinations\n");

	string output_name = "coordination.dat";
	output = (FILE *)NULL;
	output = fopen (output_name.c_str(), "w");
	if (output == (FILE *)NULL) {
		printf ("couldn't open the output file for reading!\n");
		exit (1);
	}

	printf ("\noutput file = %s\n", output_name.c_str());

	posmin = POSMIN;
	posmax = POSMAX;
	posres = POSRES;
	
	// first we figure out how many bins there are on the axis
	posbins = int ((posmax - posmin)/posres) + 1;

	axis = AXIS;

	int_low = INTERFACE_LOW;
	int_high = INTERFACE_HIGH;
	
	output_freq	=	10;					// how often the output file will be written (# of timesteps/10)
	timesteps	=	TIMESTEPS;				// # of timesteps to process through
	
	printf ("Total timesteps = %d\n", timesteps);

return;
}

void CoordinationTest::FindInterfacialWaters (std::vector<Water *>& int_mols, std::vector<Atom *>& int_atoms, AmberSystem * sys) {
	
	int_mols.clear();
	int_atoms.clear();

	Molecule * pmol;

	// go through the system
	RUN (sys->Molecules()) {
		// grab each molecule
		pmol = sys->Molecules(i);

		// we're only looking at waters for SFG analysis right now
		if (pmol->Name() != "h2o") continue;

		// first thing we do is grab a water molecule to work with
		Water * water = static_cast<Water *>(pmol);

		double position = water->Atoms(0)->Position()[axis];
		// and find molecules that sit within the interface.
		if (position < PBCFLIP) position += Atom::Size()[axis];		// adjust for funky boundaries
		// these values have to be adjusted for each system
		if (position < INTERFACE_LOW or position > INTERFACE_HIGH) continue;				// set the interface cutoffs
		
		int_mols.push_back (water);
		RUN2(water->Atoms()) {
			int_atoms.push_back (water->Atoms(j));
		}
	}

return;
}

void CoordinationTest::FindWaters (std::vector<Water *>& int_mols, std::vector<Atom *>& int_atoms, AmberSystem * sys) {
	
	int_mols.clear();
	int_atoms.clear();

	Molecule * pmol;
	Water * water;

	// go through the system
	RUN (sys->Molecules()) {
		// grab each molecule
		pmol = sys->Molecules(i);

		// we're only looking at waters for SFG analysis right now
		if (pmol->Name() != "h2o") continue;

		// first thing we do is grab a water molecule to work with
		water = static_cast<Water *>(pmol);
		// add it to the group
		int_mols.push_back (water);

		// and then add all of its atoms
		RUN2(water->Atoms()) {
			int_atoms.push_back (water->Atoms(j));
		}
	}

return;
}

int main (int argc, char **argv) {

	CoordinationTest coords;

	// set up a holder for all the different coordination types to bin them according to location in the slab - i.e. we'll make a density histogram for the different coordination types
	// The map to shove the data into the right coordination's histogram
	std::map<coordination, std::vector<int> > histograms;
	
	std::map<coordination, string>::iterator it;
	std::map<coordination, string> * mymap = &(Water::CoordinationNames);

	// here initialize all the histograms to get ready for binning the positions
	for (it = (*mymap).begin(); it != (*mymap).end(); it++) {
		
		coordination coord = it->first;
		string name = it->second;
		
		histograms[coord].resize(coords.posbins, 0);
	}

	std::vector<Water *> waters;
	std::vector<Atom *> atoms;
	for (int step = 0; step < coords.timesteps; step++) {

		// we'll go through and pick out all the waters in the system at the interface
		//coords.FindInterfacialWaters (waters, atoms, coords.sys);
		coords.FindWaters (waters, atoms, coords.sys);

		// to find the coordination of the atoms we have to update the bond graph
		coords.sys->UpdateGraph(atoms);

		// now, for each water, we find its coordination
		RUN (waters) {
			coordination coord = waters[i]->Coordination();

			// calculate its position in the slab, and find the histogram bin for it
			Atom * oxy = waters[i]->GetAtom ("O");
			VecR position = oxy->Position();
			double r = position[coords.axis];
			if (r < PBCFLIP) r += Atom::Size()[coords.axis];
			int bin = (int)((r - coords.posmin)/coords.posres);

			// then, update the respective histogram for that coordination
			histograms[coord][bin]++;
		}

		coords.sys->LoadNext();
	}


	// now that all the histograms have been compiled... let's output them to a file somehow
	vector<coordination> vc;
	vc.push_back (OH);
	vc.push_back (H);
	vc.push_back (OO);
	vc.push_back (O);
	vc.push_back (HH);
	vc.push_back (OHH);

	for (int i = 0; i < coords.posbins; i++) {
		
		double pos = double(i) * coords.posres + coords.posmin;

		printf ("% 8.3f", pos);

		RUN2 (vc) {
			coordination c = coordination(vc[j]);
			vector<int> h = histograms[c];
			double val = double(h[i]);
			printf ("% 8d", int(val/coords.timesteps));
		}

		printf ("\n");
	}

return 0;
}
