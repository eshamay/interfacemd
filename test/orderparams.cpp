#include "orderparams.h"

// output data to a file
void OrderParameters::PrintOutput () {

	rewind (output);

	if (!((step*10) % 2500)) {
		for (unsigned int pos = 0; pos < posbins; pos++) {
			double position = double(pos) * posres + posmin;

		for (unsigned int s1; s1 < angbins; s1++) {
			double S1 = double(s1)*angres + angmin;

		for (unsigned int s2n; s2n < angbins; s2n++) {
			double S2
		for (unsigned int s2d; s2d < angbins; s2d++) {

		}}}}
	}

	fflush (output);

return;
}

// a status output meter
void OrderParameters::PrintStatus () {

	if (!(timestep % 2500)) 
		cout << endl << step << ") ";
	if (!(timestep % 250))  
		cout << "*";
	
	fflush (stdout);

return;
}


int main (int argc, char **argv) {

	// start the analysis - run through each timestep
	for (timestep = 0; timestep < timesteps; timestep++) {
		// then look at each molecule
		for (unsigned int mol = 0; mol < sys->NumMols(); mol++) {

			// find all the waters
			Water * wat = static_cast<Water *>(sys->Molecules(mol));
			if (wat->Name() != "h2o") continue;

			// take each water and find its position in the slab
			double pos = wat->GetAtom("O")->Y();
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

		sys->LoadNext();

		PrintStatus();
		PrintOutput ();
	}

	PrintOutput ();

return 0;
}


