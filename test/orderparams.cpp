#include "orderparams.h"

void PrintStatus (int step);
void PrintOutput (int step, FILE * output);

int main (int argc, char **argv) {

	// here is our system for analysis
	AmberSystem sys (PRMTOP, MDCRD, FORCE);

	FILE * output = (FILE *)NULL;
	output = fopen ("OUTPUTFILE", "w");
	if (output == (FILE *)NULL) {
		printf ("couldn't open the output file for reading!\n");
		exit (1);
	}

	// here's the number of bins for the position and angle histograms
	const int posbins = (posmax-posmin)/posres;
	const int angbins = (angmax-angmin)/angres;

	// setup and initialize the system histogram
	// The histo looks like histo[y-position][S1][S2 numerator][S2 denominator]
	int histo[posbins][angbins][angbins][angbins];
	int number_density[posbins];
	for (int i=0; i<posbins; i++) {
		number_density[i] = 0;
	for (int j=0; j<angbins; j++) {
	for (int k=0; k<angbins; k++) {
	for (int l=0; l<angbins; l++) {
		histo[i][j][k][l] = 0;
	}}}


	// start the analysis - run through each timestep
	for (unsigned int step = 0; step < TIMESTEPS; step++) {
		// then look at each molecule
		for (unsigned int mol = 0; mol < sys.NumMols(); mol++) {

			// find all the waters
			Water * wat = static_cast<Water *>(sys.Molecules(i));
			if (wat->Name() != "h2o") continue;

			// take each water and find its position in the slab
			double pos = wat->GetAtom("O");
			if (pos < 15.0) pos += Atom::Size()[axis];
			int posbin = int ((pos-posmin)/posres);

			// The two order parameters are calculated from the Euler angles, so let's find those
			// first set the molecular axes up
			wat->SetOrderAxes ();
			// then find the Euler angles for the tilt and twist of the molecule
			wat->CalcEulerAngles ();
			double tilt = wat->EulerAngles[1];
			double twist = wat->EulerAngles[2];

			// calculate the S1 term
			double S1 = 3.0 * cos(tilt) * cos(tilt) - 1.0;

			// and the S2 numerator and denominator
			double S2_num = sin(tilt) * cos(2.0 * twist);
			double S2_den = sin(tilt);

			// and now bin all those values for the histogram
			int S1_bin = int ((S1-posmin)/posres);
			int S2_num_bin = int ((S2_num-angmin)/angres);
			int S2_den_bin = int ((S2_den-angmin)/angres);

			++histo[posbin][S1_bin][S2_num_bin][S2_den_bin];
			++number_density[posbin];
		}

		sys.LoadNext();

		this->PrintStatus(step);
		this->PrintOutput (step, output);

	}

return 0;
}

// output data to a file
void PrintOutput (int step, FILE * output) {

	if (!((step*10) % 2500)) {

	}

return;
}

// a status output meter
void PrintStatus (int step) {

	if (!(step % 2500)) 
		cout << endl << _timestep << ") ";
	if (!(int step % 250))  
		cout << "*";

return;
}

