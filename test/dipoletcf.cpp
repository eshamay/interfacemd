#include "dipoletcf.h"
#include "../hno3analysis.h"

void OutputStatus (int step);

int main (int *argc, char **argv) {
	
	VecR dims (xSize, ySize, zSize);
	
	XYZSystem sys (filename, dims, wannierfile);

	FILE * output = fopen (outputfilename.c_str(), "w");
	
	// run through a bunch of timesteps
	printf ("Performing a TCF analysis of the system dipole\n");
	printf ("system size for this analysis: %8.3f%8.3f%8.3f\n", dims[x], dims[y], dims[z]);
	printf ("number of timesteps = %d\n", sys.NumSteps());
	printf ("Outputting data to file: %s\n", outputfilename.c_str());

	VecR firstDipole, stepDipole;
	for (int step = 0; step < sys.NumSteps(); step++) {

		stepDipole.Zero();
int h3o = 0, no3 = 0, hno3 = 0, oh = 0, h2o = 0, und = 0;

		Molecule * mol;
		RUN (sys.Molecules()) {
			//if (sys.Molecules(i)->Name() != "h2o") continue;

			mol = sys.Molecules(i);
if (mol->Name() == "hno3") 
	hno3++;
if (mol->Name() == "h2o") 
	h2o++;
if (mol->Name() == "oh") 
	oh++;
if (mol->Name() == "no3") 
	no3++;
if (mol->Name() == "h3o")
	h3o++;
if (mol->Name() == "undefined")
	printf ("step %d\n", step);

			mol->CalcDipole();
			stepDipole += mol->Dipole();
			/*
			if (mol->Name() == "hno3" or mol->Name() == "no3") {
				printf ("% 10.3f\n", mol->Dipole().Magnitude());
			}
			*/
				//mol->Print();
		}

		if (!step)
			firstDipole = stepDipole;

		if (step)
			fprintf (output, "%f\n", stepDipole * firstDipole);

//printf ("(%d) %d %d %d %d %d %d\n", step, h2o, h3o, hno3, no3, oh);
		//OutputStatus (step);

		sys.LoadNext();
	}

	fclose(output);
		
return 0;
}

void OutputStatus (int step) {

	if (!(step % (OUTPUT_FREQ * 10))) 
		printf ("\n%d)  ", step);
	if (!(step % OUTPUT_FREQ))
		printf ("*"); fflush (stdout);

return;
}
