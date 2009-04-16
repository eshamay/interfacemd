#include "orderparams.h"

OrderParameters::OrderParameters (int argc, const char **argv, const WaterSystemParams& params)
	:	WaterSystem(argc, argv, params),
	angmax (1.0), angmin (-1.0), angres (0.01), angbins ((angmax-angmin)/angres),
	S1 (std::vector<double> (posbins, 0.0)),
	S2_num (std::vector<double> (posbins, 0.0)),
	S2_den (std::vector<double> (posbins, 0.0)),
	number_density (std::vector<int> (posbins, 0))
	
{
	
	printf ("Running an order parameter analysis with the following options:\n");

	if (argc < 3) {
		printf ("Rerun with the two interface locations:\norderparams <int_low> <int_high>\n");
		exit(1);
	}
	
	printf ("\tAngle cosines will range from:\n\t\tMax = % 8.3f\n\t\tMin = % 8.3f\n\t\tResolution = % 8.3f\n", angmin, angmax, angres);

	// The histogram looks like histo[y-position][S1][S2 numerator][S2 denominator]

return;
}
	
// output data to a file
void OrderParameters::OutputData () {

	rewind (output);

	if (!(timestep % (output_freq * 10))) {
		for (unsigned int pos = 0; pos < posbins; pos++) {

			double position = double(pos) * this->posres + this->posmin;
			double N = (double)number_density[pos];
			if (!N) continue;

			// straightforward output of the equations for the order parameters
			fprintf (output, "% 10.3f% 10.3f% 10.3f\n", 
				position,						// position along the axis
				0.5 * (3.0 * S1[pos]/N - 1), 	// S1
				S2_num[pos] / S2_den[pos]		// S2
			);
		}
	}

	fflush (output);

return;
}

/*
// a status output meter
void OrderParameters::OutputStatus () const {

	if (!(timestep % 2500)) 
		cout << endl << timestep << "/" << timesteps << " ) ";
	if (!(timestep % 250))  
		cout << "*";
	
	fflush (stdout);

return;
}
*/
void OrderParameters::Analysis () {

	// if restarting, then fast-forward to the point where we'll restart
	#ifdef RESTART
	printf ("\n\n*** Restart Run ***\n\tNow skipping %d steps before beginning analysis\n", restart);
	for (timestep = 0; timestep < restart; timestep++) {
		sys.LoadNext();
	}

	printf ("\n*** Begin Analysis ***\n\tStarting analysis at timestep %d\n\n", timestep);
	for (timestep = restart; timestep < timesteps; timestep++) {
	#else
	// start the analysis - run through each timestep
	for (timestep = 0; timestep < timesteps; timestep++) {
	#endif

		// find all the waters
		this->FindWaters ();

// ***********************************
// used for testing on so4 + no3
		this->SliceWaters (0.0, 50.0);
// ***********************************

		this->UpdateMatrix ();

		RUN (int_mols) {

			Water * wat = int_mols[i];

			// calculate its position in the slab, and find the histogram bin for it
			Atom * oxy = wat->GetAtom ("O");
			VecR r = oxy->Position();
			double position = r[axis];
			if (position < pbcflip) position += Atom::Size()[axis];		// deal with the periodic cutoffs

			// we're going to do averaging of the two interfaces.
			// First we find if the water is in the upper or lower interface and find its position relative to the gibbs dividing surface
			double distance = (position > middle) ? position - int_high : int_low - position;

			// find the position bin
			int posbin = int ((distance - this->posmin)/this->posres);

			// The two order parameters are calculated from the Euler angles, so let's find those
			// first set the molecular axes up
			wat->SetOrderAxes ();

			// then find the Euler angles for the tilt and twist of the molecule
			wat->CalcEulerAngles (this->axis);
			double tilt = wat->EulerAngles[1];
			double twist = wat->EulerAngles[2];

			// calculate the S1 term
			double S1_value = cos(tilt) * cos(tilt);

			// and the S2 numerator and denominator
			double sin_tilt_sq = sin(tilt) * sin(tilt);	// this is just done for a little more speed - not calculating sin twice... (?)
			double S2_num_value = sin_tilt_sq * cos(2.0 * twist);
			double S2_den_value = sin_tilt_sq;

/*
			// and now bin all those values for the histogram
			int S1_bin = int ((S1-posmin)/posres);
			int S2_num_bin = int ((S2_num-angmin)/angres);
			int S2_den_bin = int ((S2_den-angmin)/angres);
*/
			S1[posbin] += S1_value;
			S2_num[posbin] += S2_num_value;
			S2_den[posbin] += S2_den_value;
			++number_density[posbin];
		}

		this->sys.LoadNext();

		this->OutputStatus ();
		this->OutputData ();
	}

	this->OutputData ();

	return;
}

int main (int argc, const char **argv) {

	WaterSystemParams params;

	params.prmtop = "prmtop";
	params.mdcrd = "mdcrd";
	params.mdvel = "";
	params.axis = y;
	params.timesteps = 200000;
	params.restart = 0;
	#ifdef RESTART
		params.restart = 100000;
	#endif
	#ifdef AVG
		params.avg = true;
		params.posmin = -50.0;
		params.posmax = 50.0;
		params.output = "orderparams.avg.dat";
	#else
		params.avg = false;
		params.posmin = -5.0;
		params.posmax = 150.0;
		params.output = "orderparams.dat";
	#endif
	params.posres = 0.100;
	params.pbcflip = 15.0;
	params.output_freq = 50;

	OrderParameters par (argc, argv, params);
	
	par.Analysis();

return 0;
}
