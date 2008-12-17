#include "xyzdipole_analysis.h"

int main (int *argc, char **argv) {
	
	VecR dims (xSize, ySize, zSize);
	
	XYZSystem sys (filename, dims, wannierfile);

	FILE * output = fopen (outputfilename.c_str(), "w");
	
	// run through a bunch of timesteps
	printf ("systemsize for this analysis: %8.3f%8.3f%8.3f\n", dims[x], dims[y], dims[z]);
	printf ("number of timesteps = %d\n", sys.NumSteps());
	printf ("Outputting data to file: %s\n", outputfilename.c_str());
	for (int step = 0; step < sys.NumSteps(); step++) {
	//for (int step = 0; step < 5000; step++) {

#ifdef POWER_RUN
		double oh = 0.0, no = 0.0;
		int numH2Os = 0;

		RUN (sys.Molecules()) {
			Molecule * mol (sys.Molecules(i));

			if (mol->Name() == "h2o") {
				Water * wat = static_cast<Water *>(mol);
				wat->FindOHBonds ();
				oh += wat->OH1()->Magnitude();
				oh += wat->OH2()->Magnitude();
				++numH2Os;
			}

			if (mol->Name() == "oh") {
				Hydroxide * hyd = static_cast<Hydroxide *>(mol);
				hyd->SetAtoms();
				oh += hyd->OH()->Magnitude();
			}
			if (mol->Name() == "no3") {
				Nitrate * nit = static_cast<Nitrate *>(mol);
				nit->SetAtoms();
				no += nit->NO1()->Magnitude();
				no += nit->NO2()->Magnitude();
				no += nit->NO3()->Magnitude();
			}
		}

		fprintf (output, "%d\t%f\n", step, oh);
#endif

#ifdef TEST_RUN
//		int minID, maxID;
//		double min = 0.0, max = 0.0, dip;

		RUN (sys.Molecules()) {

			Molecule * mol (sys.Molecules(i));
			mol->CalcDipole();
			fprintf (output, "% 10.4f\n", mol->Dipole().Magnitude() * 4.8028);
		}
/*
			if (!min) {
				min = dip;
				minID = i;
				max = dip;
				maxID = i;
			}

			if (dip < min) {
				min = dip;
				minID = i;
			}

			if (dip > max) {
				max = dip;
				maxID = i;
			}
		}

		fprintf (output, "%d)\t(%s - %d) %f\t(%s - %d) %f\n", step, sys.Molecules(minID)->Name().c_str(), minID, min, sys.Molecules(maxID)->Name().c_str(), maxID, max);
*/
#endif

#ifdef TCF_RUN
		VecR firstDipole, stepDipole;
		stepDipole.Zero();

		RUN (sys.Molecules()) {
			//if (sys.Molecules(i)->Name() != "h2o") continue;
			sys.Molecules(i)->CalcDipole();
			stepDipole += sys.Molecules(i)->Dipole();
		}

		//stepDipole = sys.SystemDipole();

		if (!step)
			firstDipole = stepDipole;

		if (step)
			fprintf (output, "%f\n", stepDipole * firstDipole);

#endif

/*
		// output a status meter, and flush the output file
		int prog ((int)((double)step/numSteps*100.0));
		string progress (prog/2, '#');
		string blanks (50 - prog/2, ' ');
		ostringstream bar;
		bar << "[" << progress << prog << "%" << blanks << "]";
		printf ("%s\r", bar.str().c_str());
*/

		if (!(step % (OUTPUT_FREQ * 10))) printf ("\n%d)  ", step);
		if (!(step % OUTPUT_FREQ)) { printf ("*"); fflush (stdout); fflush (output); }
		sys.LoadNext();
	}

	fclose(output);
		
return 0;
}
