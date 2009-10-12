#include "ambertest.h"

// the command line arguments are the list of atom *NAMES* that are to be analyzed... exactly as they are named in the prmtop file
// the files prmtop, mdcrd, and mdvel should point to their respective files.
AmberTest::AmberTest (int const argc, const char **argv, const WaterSystemParams& params)
	:	WaterSystem<AmberSystem> (params)
{

	this->sys = new AmberSystem("prmtop", "mdcrd", "mdvel");

return;
}

void AmberTest::DoSomething () {

	this->FindWaters();

	Water * wat = int_wats[0];
	wat->Print();
	wat->Flip(y);
	wat->Print();

return;
}
int main (const int argc, const char **argv) {

	#ifdef AVG
		if (argc < 3) {
			printf ("for averaging, we need the two interface locations: <low> <high>\n");
			exit(1);
		}
	#endif

	WaterSystemParams params;

	params.axis = y;
	params.timesteps = 200000;
	#ifdef RESTART
		params.restart = 100000;
	#endif
	#ifdef AVG
		params.avg = true;
		params.posmin = -40.0;
		params.posmax = 40.0;
		params.output = "density.avg.100+.dat";
	#else
		params.avg = false;
		params.posmin = -20.0;
		params.posmax = 150.0;
		params.output = "density.dat";
	#endif
	params.posres = 0.100;
	params.pbcflip = 30.0;
	params.output_freq = 500;


	AmberTest test (argc, argv, params);
	test.DoSomething();

return 0;
}

