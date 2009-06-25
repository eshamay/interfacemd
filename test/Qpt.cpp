#include "Qpt.h"

int main () {

	VecR dims (xSize, ySize, zSize);

	XYZSystem sys (filename, dims, wannierfile);


	for (int step = 0; step < sys.NumSteps(); step++) {
		RUN (sys.Molecules()) {
			Molecule * na = sys.Molecules(i);
			// find the nitric acid molecule.
			if (na->Name() != "hno3") continue;

			double qpt = q_pt (sys, na);
			printf ("%10d% 10.4f\n", step, qpt);
			fflush (stdout);

		}
		sys.LoadNext();
		sys.Matrix()->UpdateMatrix();
	}


return 0;
}
