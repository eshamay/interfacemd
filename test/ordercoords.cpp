#include "ordercoords.h"

CoordOrderParams::CoordOrderParams (const int argc, const char **argv, const WaterSystemParams& params)
	:	WaterSystem(argc, argv, params),
	angmax (1.0), angmin (-1.0), angres (0.01), angbins ((angmax-angmin)/angres),
	// The histograms look like [coord_type][bin_position] = number
	S1 (Double_histo (45, Double_vec (posbins, 0.0))),
	S2_num (Double_histo (45, Double_vec (posbins, 0.0))),
	S2_den (Double_histo (45, Double_vec (posbins, 0.0))),
	number_density (Int_histo (45, Int_vec (posbins, 0.0)))

{
	
	printf ("Running an order parameter analysis on water-coordination types with the following options:\n");

	if (argc < 3) {
		printf ("Rerun with the two interface locations:\norderparams <int_low> <int_high>\n");
		exit(1);
	}
	
	printf ("\tAngle cosines will range from:\n\t\tMax = % 8.3f\n\t\tMin = % 8.3f\n\t\tResolution = % 8.3f\n", angmin, angmax, angres);

	// set up all the coordination types
	this->SetHistograms (argc, argv);

return;
}
	
void CoordOrderParams::SetHistograms (const int argc, const char **argv) {

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
 
 	// so here we grab the names of all the coordination types
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

// output data to a file
void CoordOrderParams::OutputData () {

	rewind (output);

	if (!((timestep) % (output_freq * 10))) {

		// first we output a header of all the coordination types
		fprintf (output, "%12s", "position");
		RUN (vcoords) {
			//if (vcoords[i] > HIGH_COORD) continue;
			fprintf (output, "%s\t", name_map[vcoords[i]].c_str());
		}
		fprintf (output, "\n");
	
		// The output for each row will be the position followed by S1-S2 pairs for each coordination
		for (int i = 0; i < posbins; i++) {
			
			// calculate the slab position
			double pos = double(i) * posres + posmin;

			// only print out positions that have some data in them
			int totalN = 0;			// total number density over all coordinations for this position
			RUN2 (vcoords) {
				int coord = int(vcoords[j]);
				totalN += number_density[coord][i];
			}
			if (!totalN) continue;

			// output the position for the row
			fprintf (output, "% 12.3f", pos);

			// and now go through each coordination
			RUN2 (vcoords) {
				int coord = int(vcoords[j]);

				//if (vcoords[i] > HIGH_COORD) continue;
				double N = double(number_density[coord][i]);
				//if (number_density[coord][i])
//printf ("coord = %d, i = %d, %d\n", coord, i, number_density[coord][i]);
				double s1 = 0.5 * (S1[coord][i]/N);					// the S1 order parameter
				fprintf (output, "% 10.3f", s1);

				double s2 = (S2_num[coord][i]) / (S2_den[coord][i]);		// S2 order parameter
				fprintf (output, "% 10.3f", s2);

				//printf ("% 8d% 8d% 8d% 8.3f% 8.3f% 8.3f\n", 
					//coord, int(N), i, S1[coord][i], S2_num[coord][i], S2_den[coord][i]);

			}
			fprintf (output, "\n");
		}
	}

	fflush (output);

return;
}

void CoordOrderParams::Analysis () {

	// start the analysis - run through each timestep
	for (timestep = 0; timestep < timesteps; timestep++) {

		// find all the waters in the system
		this->FindWaters ();

		this->SliceWaters (0.0, 50.0);

		// update our adjacency matrix with connectivity and bonding data
		this->UpdateMatrix ();

		RUN (int_mols) {

			Water * wat = this->int_mols[i];

			// Here we'll find out the water's coordination type (OH, OOH, OHH, etc.)
			int coord = int(this->matrix.WaterCoordination (wat));
			// note that we're only collecting data on coordinations up to a certain type...
			if (coord >= this->S1.size()) continue;
	
			// take each water and find its position in the slab
			Atom * oxy = wat->GetAtom ("O");
			VecR r = oxy->Position();
			double pos = r[axis];
			if (pos < pbcflip) pos += Atom::Size()[axis];
			// we're going to do averaging of the two interfaces.
			// First we find if the water is in the upper or lower interface and find its distance to the nearest gibbs dividing surface
			double distance = (pos > middle) ? pos - int_high : int_low - pos;

			// find the position bin
			int posbin = int ((distance-posmin)/posres);

			// we should have a data space large enough for binning, but hey, sometimes things go wacky
			if (posbin >= posbins) continue;
			
			// The two order parameters are calculated from the Euler angles, so let's find those
			// first set the molecular axes up
			wat->SetOrderAxes ();
		
			// then find the Euler angles for the tilt and twist of the molecule
			wat->CalcEulerAngles (this->axis);
			double tilt = wat->EulerAngles[1];
			double twist = wat->EulerAngles[2];
		
			// calculate the S1 term (3*cos(tilt)^2-1)
			double S1_value = 3.0 * cos(tilt) * cos(tilt) - 1.0;
		
			// and the S2 numerator and denominator (S2 = <sin(t)cos(2t)>/<sin(t)>)
			double S2_num_value = sin(tilt) * cos(2.0 * twist);
			double S2_den_value = sin(tilt);
		
			// and now to bin all that data
			this->S1[coord][posbin] += S1_value;
			this->S2_num[coord][posbin] += S2_num_value;
			this->S2_den[coord][posbin] += S2_den_value;
			this->number_density[coord][posbin] += 1;
		}
	
		this->sys.LoadNext();

		this->OutputStatus ();
		this->OutputData ();
	}

	this->OutputData ();

	return;
}

int main (const int argc, const char **argv) {

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
		params.output = "orderparams.coords.avg.dat";
	#else
		params.avg = false;
		params.posmin = -5.0;
		params.posmax = 150.0;
		params.output = "orderparams.coords.dat";
	#endif
	params.posres = 0.100;
	params.pbcflip = 15.0;
	params.output_freq = 25;

	CoordOrderParams par (argc, argv, params);
	
	par.Analysis();

return 0;
}
