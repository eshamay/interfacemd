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
	
	// here initialize all the histograms to get ready for binning the positions
	// to do this, we have to get all the different coordination types possible. Those come from the coordination type map in the bondgraph
	
	name_map = sys->bondgraph.CoordNameMap();

	coord_map::iterator coord_it, coord_end;
	vcoords.clear();
	for (coord_it = name_map.begin(); coord_it != name_map.end(); coord_it++) {

		coordination coord = coord_it->first;
		string name = coord_it->second;
		
		histo[coord].resize(posbins, 0);
		vcoords.push_back(coord);
	}

return;
}

void CoordinationTest::FindInterfacialWaters (VPWATER& int_mols, VPATOM& int_atoms) {
	
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

void CoordinationTest::FindWaters (VPWATER& int_mols, VPATOM& int_atoms) {
	
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

void CoordinationTest::OutputData (const int step) {

	if (!(step % (OUTPUT_FREQ * 10))) {

		rewind (output);

		// now that all the histograms have been compiled... let's output them to a file somehow
		// first, output a header
		fprintf (output, "position\t");
		RUN (vcoords) 
			fprintf (output, "%s\t", name_map[vcoords[i]].c_str());
		fprintf (output, "\n");
	
		// and now print all the data out
		for (int i = 0; i < posbins; i++) {
			
			double pos = double(i) * posres + posmin;
	
			fprintf (output, "% 8.3f", pos);

			RUN2 (vcoords) {
				coordination c = vcoords[j];
				int val = histo[c][i];
				fprintf (output, "% 8d", val);
			}

			fprintf (output, "\n");
		}
	}

return;
}

void CoordinationTest::OutputStatus (const int step) const {
	
	if (!(step % (OUTPUT_FREQ * 10)))
		printf ("\n%10d/%d)  ", step, TIMESTEPS);

	if (!(step % OUTPUT_FREQ)) 
		printf ("*");

	fflush (stdout);

return;
}





int main (int argc, char **argv) {

	CoordinationTest coords;

	VPWATER waters;
	VPATOM atoms;
	for (int step = 0; step < coords.timesteps; step++) {

		// we'll go through and pick out all the waters in the system at the interface
		//coords.FindInterfacialWaters (waters, atoms, coords.sys);
		coords.FindWaters (waters, atoms);

		// to find the coordination of the atoms we have to update the bond graph
		coords.sys->bondgraph.UpdateGraph(atoms);

		// now, for each water, we find its coordination
		Water * wat;
		RUN (waters) {
			wat = waters[i];
			coordination coord = coords.sys->bondgraph.WaterCoordination (wat);

			// the highest coordination that we'll look at...
			if (coord > OOHHH) continue;

			// calculate its position in the slab, and find the histogram bin for it
			Atom * oxy = wat->GetAtom ("O");
			VecR r = oxy->Position();
			double position = r[coords.axis];
			if (position < PBCFLIP) position += Atom::Size()[coords.axis];
			int bin = (int)((position - coords.posmin)/coords.posres);

			// then, update the respective histogram for that coordination
			coords.histo[coord][bin]++;
		}

		coords.sys->LoadNext();

		coords.OutputStatus (step);
		coords.OutputData (step);
	}

	coords.OutputData (0);


return 0;
}
