#include "orderparams.h"

OrderParameters::OrderParameters (int argc, char **argv) {
	
	// here is our system for analysis
	sys = new AmberSystem (PRMTOP, MDCRD, FORCE);

	printf ("Running an order parameter analysis with the following options:\n");

	string output_name = "orderparams.dat";
	output = (FILE *)NULL;
	output = fopen (output_name.c_str(), "w");
	if (output == (FILE *)NULL) {
		printf ("couldn't open the output file for reading!\n");
		exit (1);
	}

	printf ("\noutput file = %s\n", output_name.c_str());

	axis = y;
	
	output_freq	=	10;					// how often the output file will be written (# of timesteps/10)
	timesteps	=	200000;				// # of timesteps to process through
	
	printf ("timesteps = %d\n", timesteps);

	// position boundaries and bin width
	posmin	= -25.0;
	posmax	= 25.0;
	posres	= 0.5;
	posbins = (posmax-posmin)/posres;

	printf ("Position analysis will range from % 8.3f to % 8.3f angstroms of the gibbs dividing surface, in % 8.3f increments\n", posmin, posmax, posres);

	if (argc < 3) {
		printf ("Rerun with the two interface locations:\norderparams <int_low> <int_high>\n");
		exit(1);
	}
	int_low = atof(argv[1]);
	int_high = atof(argv[2]);
	middle = (int_low + int_high)/2.0;

	printf ("Low interface : % 8.3f\nHigh interface : % 8.3f\n", int_low, int_high);
	
	angmax	 = 1.0;
	angmin	 = -1.0;
	angres	 = 0.01;
	angbins = (angmax-angmin)/angres;
	
	printf ("Angle cosines will be measured in increments of % 8.3f\n\n", angres);

	// setup and initialize the system histogram
	// The histo looks like histo[y-position][S1][S2 numerator][S2 denominator]
	S1.clear(); S1.resize (posbins, 0.0);
	S2_num.clear(); S2_num.resize (posbins, 0.0);
	S2_den.clear(); S2_den.resize (posbins, 0.0);
	number_density.clear(); number_density.resize(posbins, 0);

return;
}
	
// output data to a file
void OrderParameters::PrintOutput () {

	rewind (output);

	if (!((timestep*10) % 2500)) {
		for (unsigned int pos = 0; pos < posbins; pos++) {

			double position = double(pos) * posres + posmin;
			double N = (double)number_density[pos];
			if (N == 0) continue;

			fprintf (output, "% 10.3f% 10.3f% 10.3f\n", 
				position,
				0.5 * S1[pos]/N,
				S2_num[pos] / S2_den[pos]
			);
		}
	}

	fflush (output);

return;
}

// a status output meter
void OrderParameters::PrintStatus () {

	if (!(timestep % 2500)) 
		cout << endl << timestep << ") ";
	if (!(timestep % 250))  
		cout << "*";
	
	fflush (stdout);

return;
}


int main (int argc, char **argv) {

	OrderParameters par (argc, argv);

	// start the analysis - run through each timestep
	for (par.timestep = 0; par.timestep < par.timesteps; par.timestep++) {
		// then look at each molecule
		for (unsigned int mol = 0; mol < par.sys->NumMols(); mol++) {

			// find all the waters
			Water * wat = static_cast<Water *>(par.sys->Molecules(mol));
			if (wat->Name() != "h2o") continue;

			// take each water and find its position in the slab
			double pos = wat->GetAtom("O")->Y();
			if (pos < 15.0) pos += Atom::Size()[par.axis];

			// we're going to do averaging of the two interfaces.
			// First we find if the water is in the upper or lower interface and find its position relative to the gibbs dividing surface
			double distance = (pos > par.middle) ? pos - par.int_high : par.int_low - pos;

			int posbin = int ((distance-par.posmin)/par.posres);

			// The two order parameters are calculated from the Euler angles, so let's find those
			// first set the molecular axes up
			wat->SetOrderAxes ();
			// then find the Euler angles for the tilt and twist of the molecule
			wat->CalcEulerAngles (par.axis);
			double tilt = wat->EulerAngles[1];
			double twist = wat->EulerAngles[2];

			// calculate the S1 term
			double S1_value = 3.0 * cos(tilt) * cos(tilt) - 1.0;

			// and the S2 numerator and denominator
			double S2_num_value = sin(tilt) * cos(2.0 * twist);
			double S2_den_value = sin(tilt);

/*
			// and now bin all those values for the histogram
			int S1_bin = int ((S1-posmin)/posres);
			int S2_num_bin = int ((S2_num-angmin)/angres);
			int S2_den_bin = int ((S2_den-angmin)/angres);
*/

			par.S1[posbin] += S1_value;
			par.S2_num[posbin] += S2_num_value;
			par.S2_den[posbin] += S2_den_value;
			++par.number_density[posbin];
		}

		par.sys->LoadNext();

		par.PrintStatus ();
		par.PrintOutput ();
	}

	par.PrintOutput ();

return 0;
}
