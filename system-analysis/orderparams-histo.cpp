#include "orderparams-histo.h"

OrderParameters::OrderParameters (int argc, const char **argv, const WaterSystemParams& params)
	:	WaterSystem<AmberSystem>(params),
	angmax (1.0), angmin (-1.0), angres (0.01), angbins ((angmax-angmin)/angres),
	_data (std::vector< std::vector<int> > (angbins, std::vector<int> (posbins, 0)))
{

	printf ("Running an order parameter analysis to generate a histogram data matrix with the following options:\n");

	#ifdef AVG
	if (argc < 3) {
		printf ("Rerun with the two interface locations:\norderparams <int_low> <int_high>\n");
		exit(1);
	}
	#endif

	printf ("\tAngle cosines will range from:\n\t\tMax = % 8.3f\n\t\tMin = % 8.3f\n\t\tResolution = % 8.3f\n", angmin, angmax, angres);

	this->sys = new AmberSystem("prmtop", "mdcrd", "mdvel");

	// some info before starting
	printf ("The AMBER system loaded contains %d atoms\n", sys->size());

	// The histogram looks like histo[y-position][S1][S2 numerator][S2 denominator]

return;
}

// output data to a file
void OrderParameters::OutputData () {

	rewind (output);

	if (!(timestep % (output_freq * 10))) {
		// angle value is the row
		for (int ang = 0; ang < angbins; ang++) {

			// position is the column
			for (int pos = 0; pos < posbins; pos++) {

				// print out the position for each position column
				fprintf (output, "% 12d", _data[ang][pos]);

			// and then advance to the next row
			}
			fprintf (output, "\n");
		}
	}

	fflush (output);

return;
}

void OrderParameters::Analysis () {

	// if restarting, then fast-forward to the point where we'll restart
	#ifdef RESTART
	printf ("\n\n*** Restart Run ***\n\tNow skipping %d steps before beginning analysis\n", restart);
	for (timestep = 0; timestep < restart; timestep++) {
		sys->LoadNext();
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
		//this->SliceWaters (0.0, 50.0);
// ***********************************

		/******
		 * When doing any work involving H-bonding or bond distances...
		 ******/
		//this->UpdateMatrix ();
		VecR ref_axis (0.0, 1.0, 0.0);	// reference axis normal to the interfaces

		RUN (int_wats) {

			Water * wat = int_wats[i];

			// calculate its position in the slab, and find the histogram bin for it
			Atom * oxy = wat->GetAtom ("O");
			VecR r = oxy->Position();
			double position = r[axis];
			if (position < pbcflip) position += Atom::Size()[axis];		// deal with the periodic cutoffs

		#ifdef AVG
			// we're going to do averaging of the two interfaces.
			// First we find if the water is in the upper or lower interface and find its position relative to the gibbs dividing surface
			double distance = (position > middle) ? position - int_high : int_low - position;

			// find the position bin
			int posbin = int ((distance - this->posmin)/this->posres);
		#else
			int posbin = int ((position - this->posmin)/this->posres);
		#endif


			// The two order parameters are calculated from the Euler angles, so let's find those
			// first set the molecular axes up
			wat->SetOrderAxes ();

			// calculate the angles theta_z and theta_y - between the respective molecular axes and the system reference axis
			//double cos_theta_z = (wat->Z() < ref_axis);
			double cos_theta_y = (wat->Y() < ref_axis);

			// the perpendicular angle has to be treated, however, because the vector may only have values between 0 and 90 degrees... thus:
			//double theta_y = acos(cos_theta_y);
			cos_theta_y = fmod (cos_theta_y, 1.0);
			//theta_y = fmod (theta_y, M_PI/2.0);
			//cos_theta_y = cos(theta_y);
			//double bisector = (wat->Bisector() < ref_axis);

/*
			VecR oh1 = *wat->OH1();
			VecR oh2 = *wat->OH2();

			double cos_oh1 = (oh1 < ref_axis);
			double cos_oh2 = (oh2 < ref_axis);

			// The oh-vectors are not considered to be equivalent... the one with the smaller angle - more perpendicular to the surface normal - is OH-A.
			double cos_a = (cos_oh1 < cos_oh2) ? cos_oh1 : cos_oh2;
			double cos_b = (cos_oh1 > cos_oh2) ? cos_oh1 : cos_oh2;
*/
			double angbin = int((cos_theta_y - this->angmin)/this->angres);
			_data[angbin][posbin]++;
		}

		this->sys->LoadNext();

		this->OutputStatus ();
		this->OutputData ();
	}

	this->OutputData ();

	return;
}

int main (int argc, const char **argv) {

	WaterSystemParams params;

	params.axis = y;
	params.timesteps = 200000;
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
		params.posmin = -20.0;
		params.posmax = 150.0;
		params.output = "orderparams.normal.histo.dat";
	#endif
	params.posres = 0.100;
	params.pbcflip = 30.0;
	params.output_freq = 500;

	OrderParameters par (argc, argv, params);

	par.Analysis();

return 0;
}
