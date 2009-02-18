#include "coordination.h"

CoordinationTest::CoordinationTest () : 
	sys(AmberSystem(PRMTOP, MDCRD, FORCE)),
	output((FILE *)NULL),
	posmin(POSMIN), posmax(POSMAX), posres(POSRES),
	posbins(int ((posmax - posmin)/posres) + 1),
	axis(AXIS),
	int_low(INTERFACE_LOW), int_high(INTERFACE_HIGH),
	output_freq(10),
	timesteps(TIMESTEPS),
	histo (COORD_HISTOGRAM (44, HISTOGRAM (posbins, 0)))
{
	
	// here is our system for analysis
	//sys = new AmberSystem (PRMTOP, MDCRD, FORCE);

	printf ("Running a test to find water coordinations\n");

	string output_name = "coordination.dat";
	output = fopen (output_name.c_str(), "w");
	if (output == (FILE *)NULL) {
		printf ("couldn't open the output file for reading!\n");
		exit (1);
	}
	printf ("\noutput file = %s\n", output_name.c_str());

	printf ("Total timesteps = %d\n", timesteps);
	
	// here initialize all the histograms to get ready for binning the positions
	// to do this, we have to get all the different coordination types possible. Those come from the coordination type map in the bondgraph
	
 	name_map[UNBOUND] = "UNBOUND";
 	name_map[O] = "O";
 	name_map[OO] = "OO";
 	name_map[OOO] = "OOO";
 	name_map[H] = "H";
 	name_map[OH] = "OH";
 	name_map[OOH] = "OOH";
 	name_map[OOOH] = "OOOH";
 	name_map[HH] = "HH";
 	name_map[OHH] = "OHH";
 	name_map[OOHH] = "OOHH";
 	name_map[OOOHH] = "OOOHH";
 	name_map[HHH] = "HHH";
 	name_map[OHHH] = "OHHH";
 	name_map[OOHHH] = "OOHHH";
 	name_map[OOOHHH] = "OOOHHH";
 

	coord_map::iterator coord_it, coord_end;
	vcoords.clear();
	for (coord_it = name_map.begin(); coord_it != name_map.end(); coord_it++) {

		coordination coord = coord_it->first;
		string name = coord_it->second;
		
		//histo[coord].resize(posbins, 0);
		vcoords.push_back(coord);
	}

return;
}

void CoordinationTest::FindInterfacialWaters (Water_ptr_vec& int_mols, Atom_ptr_vec& int_atoms) {
	
	int_mols.clear();
	int_atoms.clear();

	Molecule * pmol;

	// go through the system
	for (int i = 0; i < sys.NumMols(); i++) {
		// grab each molecule
		pmol = sys.Molecules(i);

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

void CoordinationTest::FindWaters (Water_ptr_vec& int_mols, Atom_ptr_vec& int_atoms) {
	
	int_mols.clear();
	int_atoms.clear();

	Molecule * pmol;
	Water * water;

	// go through the system
	for (int i = 0; i < sys.NumMols(); i++) {
		// grab each molecule
		pmol = sys.Molecules(i);

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
		RUN (vcoords) {
			//if (vcoords[i] > HIGH_COORD) continue;
			fprintf (output, "%s\t", name_map[vcoords[i]].c_str());
		}
		fprintf (output, "\n");
	
		// and now print all the data out
		for (int i = 0; i < posbins; i++) {
			
			double pos = double(i) * posres + posmin;
	
			fprintf (output, "% 8.3f", pos);

			RUN2 (vcoords) {
				coordination c = vcoords[j];
				//if (vcoords[i] > HIGH_COORD) continue;
				int val = histo[(int)c][i];
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

	Water_ptr_vec waters;
	Atom_ptr_vec atoms;
	for (int step = 0; step < coords.timesteps; step++) {

		// we'll go through and pick out all the waters in the system at the interface
		//coords.FindInterfacialWaters (waters, atoms, coords.sys);
		coords.FindWaters (waters, atoms);

		// to find the coordination of the atoms we have to update the bond graph
		coords.matrix.UpdateMatrix(atoms);

		// now, for each water, we find its coordination
		Water * wat;
		RUN (waters) {
			wat = waters[i];
			coordination coord = coords.matrix.WaterCoordination (wat);

			// the highest coordination that we'll look at...
			//if (coord > HIGH_COORD) continue;

			// calculate its position in the slab, and find the histogram bin for it
			Atom * oxy = wat->GetAtom ("O");
			VecR r = oxy->Position();
			double position = r[coords.axis];
			if (position < PBCFLIP) position += Atom::Size()[coords.axis];
			int bin = (int)((position - coords.posmin)/coords.posres);

			// then, update the respective histogram for that coordination
			coords.histo[(int)coord][bin]++;
		}

		coords.sys.LoadNext();

		coords.OutputStatus (step);
		coords.OutputData (step);
	}

	coords.OutputData (0);


return 0;
}
