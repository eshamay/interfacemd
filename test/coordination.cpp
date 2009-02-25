#include "coordination.h"

CoordinationTest::CoordinationTest (int argc, char **argv, const WaterSystemParams& params) : 
	WaterSystem(argc, argv, params),
	histo (COORD_HISTOGRAM (45, HISTOGRAM (posbins, 0)))
{
	
	// here is our system for analysis
	//sys = new AmberSystem (PRMTOP, MDCRD, FORCE);

	printf ("***Data Analysis***\nRunning a test to find water coordinations\n");

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

void CoordinationTest::OutputStatus () const {
	
	if (!(timestep % (output_freq * 10)))
		printf ("\n%10d/%d)  ", timestep, timesteps);

	if (!(timestep % output_freq)) 
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
	for (timestep = 0; timestep < RESTART; timestep++) {
		sys.LoadNext();
	}

	for (timestep = RESTART; timestep < timesteps; timestep++) {
	#else 
	for (timestep = 0; timestep < timesteps; timestep++) {
	#endif

		// we'll go through and pick out all the waters in the system at the interface
		//this->FindInterfacialWaters (waters, atoms, coords.sys);
		this->FindWaters ();

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

int main (int argc, char **argv) {

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
	params.restart = 0;
	#ifdef AVG
		params.avg = true;
		params.posmin = -40.0;
		params.posmax = 40.0;
		params.output = "coordination.avg.dat";
	#else
		params.avg = false;
		params.posmin = -5.0;
		params.posmax = 150.0;
		params.output = "coordination.dat";
	#endif
	params.posres = 0.100;
	params.pbcflip = 15.0;
	params.output_freq = 25;

	CoordinationTest coords (argc, argv, params);

	coords.Analysis ();

return 0;
}
