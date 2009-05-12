#include "coordination.h"

CoordinationTest::CoordinationTest (const int argc, const char **argv, const WaterSystemParams& params) :
	WaterSystem(argc, argv, params),
	histo (Int_histo (45, Int_vec (posbins, 0)))
{

	// here is our system for analysis

	printf ("***Data Analysis***\nRunning a test to find water coordinations\n");

	this->InitCoordMaps();

return;
}

void CoordinationTest::InitCoordMaps () {

	// get all the different coordination types possible. Those come from the coordination type map in the bondgraph

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

void CoordinationTest::OutputData () {

	if (!(timestep % (output_freq * 10))) {

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

void CoordinationTest::BinPosition (Water const * const wat, coordination const coord) {

	// the highest coordination that we'll look at...
	if (coord < HIGH_COORD) {

		// calculate its position in the slab, and find the histogram bin for it
		Atom * oxy = wat->GetAtom ("O");
		VecR r = oxy->Position();
		double position = r[axis];

		if (position < pbcflip) position += Atom::Size()[axis];		// deal with the periodic cutoffs

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

void CoordinationTest::Analysis () {

	// if restarting, then fast-forward to the point where we'll restart
	#ifdef RESTART
	printf ("\n\n*** Restart Run ***\n\tNow skipping %d steps before beginning analysis\n", restart);
	for (timestep = 0; timestep < restart; timestep++) {
		sys.LoadNext();
	}

	printf ("\n*** Begin Analysis ***\n\tStarting analysis at timestep %d\n\n", timestep);
	for (timestep = restart; timestep < timesteps; timestep++) {
	#else
	for (timestep = 0; timestep < timesteps; timestep++) {
	#endif

		// we'll go through and pick out all the waters in the system at the interface
		//this->FindInterfacialWaters (waters, atoms, coords.sys);
		this->FindWaters ();

// For testing on Na2SO4 and maybe NaNO3
		this->SliceWaters (0.0, 50.0);

		// to find the coordination of the atoms we have to update the bond graph
		this->UpdateMatrix();

		// now, for each water, we find its coordination
		Water * wat;
		RUN (int_mols) {
			wat = int_mols[i];

			coordination coord = this->matrix.WaterCoordination (wat);
			this->BinPosition (wat, coord);
		}

		this->sys.LoadNext();

		this->OutputStatus ();
		this->OutputData ();
	}

	//final output
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
		params.restart = 100000;
	#endif
	#ifdef AVG
		params.avg = true;
		params.posmin = -40.0;
		params.posmax = 40.0;
		params.output = "coordination.avg.100+.dat";
	#else
		params.avg = false;
		params.posmin = -5.0;
		params.posmax = 150.0;
		params.output = "coordination.dat";
	#endif
	params.posres = 0.100;
	params.pbcflip = 15.0;
	params.output_freq = 100;

	CoordinationTest coords (argc, argv, params);

	coords.Analysis ();

return 0;
}
