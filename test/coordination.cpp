#include "coordination.h"

CoordinationTest::CoordinationTest (int argc, char **argv) : 
	sys(AmberSystem(PRMTOP, MDCRD, FORCE)),
	output((FILE *)NULL),
	posmin(POSMIN), posmax(POSMAX), posres(POSRES), middle((int_low + int_high)/2.0),
	posbins(int ((posmax - posmin)/posres) + 1),
	axis(AXIS),
	#ifdef AVG
	int_low(atof(argv[1])), int_high(atof(argv[2])),
	#endif
	output_freq(10),
	timesteps(TIMESTEPS),
	histo (COORD_HISTOGRAM (45, HISTOGRAM (posbins, 0)))
{
	
	// here is our system for analysis
	//sys = new AmberSystem (PRMTOP, MDCRD, FORCE);

	printf ("Running a test to find water coordinations\n");

#ifdef AVG
	printf ("Averaging the two interfaces - make sure the interfaces are located as follows:\n***  low interface = % 5.3f\n***  high interface = % 5.3f\n", int_low, int_high);
	string output_name = "coordination.avg.dat";
#else
	string output_name = "coordination.dat";
#endif

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
		if (position < posmin or position > posmax) continue;				// set the interface cutoffs
		
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

void CoordinationTest::BinPosition (Water const * const wat, coordination const coord) {

	// the highest coordination that we'll look at...
	if (coord < HIGH_COORD) {

		// calculate its position in the slab, and find the histogram bin for it
		Atom * oxy = wat->GetAtom ("O");
		VecR r = oxy->Position();
		double position = r[axis];

		if (position < PBCFLIP) position += Atom::Size()[axis];		// deal with the periodic cutoffs

	#ifdef AVG
		double distance = (position > middle) ? position - int_high : int_low - position;
		int bin = (int)((distance - posmin)/posres);
	#else
		int bin = (int)((position - posmin)/posres);
	#endif

		// then, update the respective histogram for that coordination
		//std::cout << bin << std::endl;
		histo[(int)coord][bin]++;
	}
return;
}



int main (int argc, char **argv) {

	#ifdef AVG
		if (argc < 3) {
			printf ("for averaging, we need the two interface locations: <low> <high>\n");
			exit(1);
		}
	#endif
	CoordinationTest coords (argc, argv);

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
			coords.BinPosition (wat, coord);
		}

		coords.sys.LoadNext();

		coords.OutputStatus (step);
		coords.OutputData (step);
	}

	coords.OutputData (0);


return 0;
}
