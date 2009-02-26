#include "watersystem.h"

WaterSystem::WaterSystem (const WaterSystemParams& params) :
	sys(AmberSystem(params.prmtop, params.mdcrd, params.mdvel)),
	output(fopen(params.output.c_str(), "w")),
	posmin(params.posmin), posmax(params.posmax), posres(params.posres), pbcflip(params.pbcflip),
	posbins(int ((posmax - posmin)/posres) + 1),
	axis(params.axis),
	output_freq(params.output_freq),
	timesteps(params.timesteps), restart(params.restart)
{
	
	if (params.avg) {	// when averaging, the interface locations are taken from the command line.
		printf ("WaterSystem::WaterSystem () --> Interface locations are needed\n");
		exit(1);
	}
	else {
		int_low = 0.0;
		int_high = 0.0;
		middle = 0.0;
	}

	OpenFile ();
	OutputHeader(params);

	return;
}

WaterSystem::WaterSystem (const int argc, const char **argv, const WaterSystemParams& params) : 
	sys(AmberSystem(params.prmtop, params.mdcrd, params.mdvel)),
	output(fopen(params.output.c_str(), "w")),
	posmin(params.posmin), posmax(params.posmax), posres(params.posres), pbcflip(params.pbcflip),
	posbins(int ((posmax - posmin)/posres) + 1),
	axis(params.axis),
	output_freq(params.output_freq),
	timesteps(params.timesteps), restart(params.restart)
	
{
	if (params.avg) {	// when averaging, the interface locations are taken from the command line.
		int_low = atof(argv[1]);
		int_high = atof(argv[2]);
		middle = (int_low + int_high)/2.0;
	}
	else {
		int_low = 0.0;
		int_high = 0.0;
		middle = 0.0;
	}

	OpenFile ();
	OutputHeader(params);

	return;
}

WaterSystem::~WaterSystem () {
	fclose (output);

	return;
}

void WaterSystem::OpenFile () {
	
	//output = fopen (params, 'w');
	if (output == (FILE *)NULL) {
		printf ("WaterSystem::WaterSystem (argc, argv) - couldn't open the output file!\n");
		exit(1);
	}

	return;
}

void WaterSystem::OutputHeader(const WaterSystemParams& params) const {

	printf ("Analysis Parameters:\n\tprmtop = \"%s\"\n\tmdcrd  = \"%s\"\n\tmdvel  = \"%s\"\n\n\tOutput Filename = \"%s\"\n\tScreen output frequency = 1/%d\n\n\tPosition extents for analysis:\n\t\tMin = % 8.3f\n\t\tMax = % 8.3f\n\t\tPosition Resolution = % 8.3f\n\n\tPrimary Axis = %d\nNumber of timesteps to be analyzed = %d\n",
			params.prmtop.c_str(), params.mdcrd.c_str(), params.mdvel.c_str(), params.output.c_str(), params.output_freq, params.posmin, params.posmax, params.posres, int(params.axis), params.timesteps);

	if (params.avg) {
		printf ("\n\nThe analysis is averaging about the two interfaces located as:\n\tLow  = % 8.3f\n\tHigh = % 8.3f\n\n", int_low, int_high);
	}

	return;
}

void WaterSystem::FindInterfacialWaters () {
	
	int_mols.clear();
	int_atoms.clear();

	Molecule * pmol;

	// go through the system
	for (int i = 0; i < sys.NumMols(); i++) {
		// grab each molecule
		pmol = sys.Molecules(i);

		// we're only looking at waters for SFG analysis right now
		if (pmol->Name() != "h2o") continue;

		// first thing we do is grab a water molecule to work with
		Water * water = static_cast<Water *>(pmol);

		double position = water->Atoms(0)->Position()[axis];
		// and find molecules that sit within the interface.
		if (position < pbcflip) position += Atom::Size()[axis];		// adjust for funky boundaries
		// these values have to be adjusted for each system
		if (position < posmin or position > posmax) continue;				// set the interface cutoffs
		
		int_mols.push_back (water);
		RUN2(water->Atoms()) {
			int_atoms.push_back (water->Atoms(j));
		}
	}

return;
}

void WaterSystem::FindWaters () {
	
	int_mols.clear();
	int_atoms.clear();

	Molecule * pmol;
	Water * water;

	// go through the system
	for (int i = 0; i < sys.NumMols(); i++) {
		// grab each molecule
		pmol = sys.Molecules(i);

		// we're only looking at waters for SFG analysis right now
		if (pmol->Name() != "h2o") continue;

		// first thing we do is grab a water molecule to work with
		water = static_cast<Water *>(pmol);
		// add it to the group
		int_mols.push_back (water);

		// and then add all of its atoms
		RUN2(water->Atoms()) {
			int_atoms.push_back (water->Atoms(j));
		}
	}

return;
}

void WaterSystem::OutputStatus () const {

	if (!(timestep % (output_freq * 10))) 
		cout << endl << timestep << "/" << timesteps << " ) ";
	if (!(timestep % output_freq))  
		cout << "*";
	
	fflush (stdout);

return;
}
