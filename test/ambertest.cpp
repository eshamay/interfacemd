#include "ambertest.h"

// the command line arguments are the list of atom *NAMES* that are to be analyzed... exactly as they are named in the prmtop file
// the files prmtop, mdcrd, and mdvel should point to their respective files.
AmberTest::AmberTest (int const argc, const char **argv, const WaterSystemParams& params)
	:	AmberSystem ("prmtop", "mdcrd", "mdvel")
{

	for (int step = 0; step < params.timesteps; step++) {
		printf ("time = %d\n", step);
		RUN (_atoms) {
			Atom * pa = this->Atoms(0);
			if (pa->Name() != "O") continue;
			if (pa->Position()[y] == 0.000)
				pa->Print();
		}
		this->LoadNext();
	}

return;
}

int main (const int argc, const char **argv) {

	WaterSystemParams params;

	params.axis = y;
	params.timesteps = 40000;
	params.avg = false;
	params.posmin = -20.0;
	params.posmax = 150.0;
	params.posres = 0.100;
	params.pbcflip = 15.0;
	params.output_freq = 100;

	AmberTest test (argc, argv, params);


return 0;
}

