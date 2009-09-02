#include "orderparams.h"

OrderParameters::OrderParameters (int argc, const char **argv, const WaterSystemParams& params)
	:	WaterSystem<AmberSystem>(params),
	angmax (1.0), angmin (-1.0), angres (0.01), angbins ((angmax-angmin)/angres),
	_data (std::vector< std::vector<double> > (25, std::vector<double> (posbins, 0.0))),
	number_density (std::vector<unsigned long int> (posbins, 0))
	/*
	phi (std::vector<double> (posbins, 0.0)),
	theta (std::vector<double> (posbins, 0.0)),
	psi (std::vector<double> (posbins, 0.0)),
	S1 (std::vector<double> (posbins, 0.0)),
	S2_num (std::vector<double> (posbins, 0.0)),
	S2_den (std::vector<double> (posbins, 0.0)),
	*/

{

	printf ("Running an order parameter analysis with the following options:\n");

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
		for (int pos = 0; pos < posbins; pos++) {

			double position = double(pos) * this->posres + this->posmin;
			unsigned long int N = number_density[pos];
			if (!N) continue;

			/*
			// straightforward output of the equations for the order parameters
			fprintf (output, "% 10.3f% 10.3f% 10.3f\n",
				position,						// position along the axis
				0.5 * (3.0 * S1[pos]/N - 1.0), 	// S1
				S2_num[pos] / S2_den[pos]		// S2
			);
			*/

			// print out the position for each data row
			fprintf (output, "% 12.5f", position);

			// then each bit of compiled data
			for (int i = 0; i < (int)_data.size(); i++) {
				fprintf (output, "% 15.5f", _data[i][pos]);
			}

			// and lastly the number density for that location in the slab
			fprintf (output, "% 15d\n", (int)N);
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

		RUN (int_mols) {

			Water * wat = int_mols[i];

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

			// then find the Euler angles for the molecule
			wat->CalcEulerAngles (this->axis);

			double phi_val = wat->EulerAngles[0];
			double theta_val = wat->EulerAngles[1];
			// for theta, if the molecule is on the bottom interface we need to adjust for the fact that the normal points in the opposite direction of the top interface. This is done just by adding 180 degrees (pi) to the angle value
			//if (position < middle)
				//theta_val += M_PI;
			double psi_val = wat->EulerAngles[2];

			// calculate the S1 term
			//double S1_value = cos(theta_val) * cos(theta_val);

			// and the S2 numerator and denominator
			//double S2_num_value = sin(theta_val) * sin(theta_val) * cos(2.0 * psi_val);
			//double S2_den_value = sin(theta_val) * sin(theta_val);

			// bin the euler angle values (cos() of...)
			_data[0][posbin] += phi_val;							// direct angle values
			_data[1][posbin] += theta_val;
			_data[2][posbin] += psi_val;
			_data[3][posbin] += sin(phi_val);						// <sin()>
			_data[4][posbin] += sin(theta_val);
			_data[5][posbin] += sin(psi_val);
			_data[6][posbin] += sin(2.0*phi_val);					// <sin(2*angle)>
			_data[7][posbin] += sin(2.0*theta_val);
			_data[8][posbin] += sin(2.0*psi_val);
			_data[9][posbin] += sin(phi_val)*sin(phi_val);			// <sin^2()>
			_data[10][posbin] += sin(theta_val)*sin(theta_val);
			_data[11][posbin] += sin(psi_val)*sin(psi_val);
			_data[12][posbin] += cos(phi_val);						// <cos()>
			_data[13][posbin] += cos(theta_val);
			_data[14][posbin] += cos(psi_val);
			_data[15][posbin] += cos(2.0*phi_val);					// <cos(2.0*angle)>
			_data[16][posbin] += cos(2.0*theta_val);
			_data[17][posbin] += cos(2.0*psi_val);
			_data[18][posbin] += cos(phi_val)*cos(phi_val);			// <cos^2()>
			_data[19][posbin] += cos(theta_val)*cos(theta_val);
			_data[20][posbin] += cos(psi_val)*cos(psi_val);
			_data[21][posbin] += sin(theta_val)*cos(2.0*phi_val);					// <sin(theta)*cos(2phi)>
			_data[22][posbin] += sin(theta_val)*cos(2.0*psi_val);					// <sin(theta)*cos(2psi)>
			_data[23][posbin] += sin(theta_val)*sin(theta_val)*cos(2.0*phi_val);	// <sin^2(theta)*cos(2phi)>
			_data[24][posbin] += sin(theta_val)*sin(theta_val)*cos(2.0*psi_val);	// <sin^2(theta)*cos(2psi)>
			++number_density[posbin];
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
		params.posmin = -5.0;
		params.posmax = 150.0;
		params.output = "orderparams.dat";
	#endif
	params.posres = 0.100;
	params.pbcflip = 15.0;
	params.output_freq = 100;

	OrderParameters par (argc, argv, params);

	par.Analysis();

return 0;
}
