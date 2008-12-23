#include "orderparams.h"

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

	OrderParameters par;

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
			int posbin = int ((pos-par.posmin)/par.posres);

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
