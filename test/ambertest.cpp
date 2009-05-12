#include <stdlib.h>
#include <sstream>
#include "ambertest.h"

int main (int argc, char *argv) {

	AmberSystem sys (PRMTOP, MDCRD, FORCE);

	const int steps = 10000;

	for (int step = 0; step < steps; step++) {

		if (!(steps % 1000)) {
		int prog = (int)((double)step / (double)steps * 100.0);
		string progress (prog/2, '#');
		string blanks (50 - prog/2, ' ');
		ostringstream percent;
		percent << prog;
		string bar = "[" + progress + percent.str() + "%" + blanks + "]";
		printf ("%s\r", bar.c_str());
		}
	}


/*
	Molecule * mol = sys.Molecules(0);

	Atom * o = (*mol)["O"];
	Atom * h = (*mol)["H2"];
	//printf ("charge) O=%10.4f\tH=%10.4f\n", o->Charge(), h->Charge());

	Molecule oh;
	oh.AddAtom (o);
	oh.AddAtom (h);
	oh.Print();
	printf ("%f\n", (o->Position() - h->Position()).Magnitude());
*/

/*
	RUN (mol->Atoms()) {
		Atom * atom = mol->Atoms()[i];

		RUN2 (atom->HBonds()) {
			Atom * hbond = atom->HBonds()[j];

			printf ("%5d(%s) : %5d(%s)\t%10.4f\n",
				atom->ID(), atom->Name().c_str(), hbond->ID(), hbond->Name().c_str(), *atom - *hbond);
		}
	}
*/

return 0;
}
