#include "densitytest.h"

// the command line arguments are the list of atom *NAMES* that are to be analyzed... exactly as they are named in the prmtop file
// the files prmtop, mdcrd, and mdvel should point to their respective files.
DensityAnalyzer::DensityAnalyzer (char * argv[], int const argc, int const numSteps, double const start, double const end, double const binsize, coord const axis) : _sys(AmberSystem(PRMTOP, MDCRD, FORCE)) {


	// this will be our output data file
	_output = (FILE *)NULL;
	string outputfile;
	#ifdef AVG
	outputfile = "density.avg.dat";
	#else
	outputfile = "density.dat";
	#endif
	_output = fopen (outputfile.c_str(), "w");
	if (_output == (FILE *)NULL) {
		printf ("couldn't open the output file density.dat. Now Exiting\n");
		exit (1);
	}
	printf ("Performing a DENSITY analysis of the system\n");
	printf ("Data will be output to the file: %s\n", outputfile.c_str());

	_step		= 0;
	_steps 		= numSteps;
	_start 		= start;
	_end 		= end;
	_binsize 	= binsize;
	_axis		= axis;

	// first we figure out how many bins there are on each axis
	_posbins = int ((_end - _start)/_binsize) + 1;

	// the location of the two interfaces
	int_low = INT_LOW;
	int_high = INT_HIGH;
	middle = (int_low + int_high)/2.0;
	#ifdef AVG
	printf ("Averaging the two interfaces for the final output\n");
	printf ("Interfaces are set at % 8.3f and % 8.3f **double-check this**\n", int_low, int_high);
	#else
	printf ("Not averaging the two interfaces - output will show total system\n");
	#endif
	printf ("Position bounds are set at: % 8.3f to % 8.3f in % 8.3f increments\n", _start, _end, _binsize);

	// from the command line, grab all the atom name that we'll be working with
	for (int i = 1; i < argc; i++) {
		_atomNames.push_back (argv[i]);
	}
	
	// some info before starting
	printf ("Analyzing densities of the atoms: ");
	RUN (_atomNames) {
		printf ("%s  ", _atomNames[i].c_str());
	}
	printf ("\n");

	// now we set up the histogram(s) such that we may bin the positions of each atom
	_density.resize(_atomNames.size(), vector<int> (_posbins, 0));
	
#ifdef DEBUG
this->Debug ("Constructor\n"); 
#endif

return;
}

void DensityAnalyzer::_PrintToFile () {

#ifdef DEBUG
this->Debug ("DensityAnalyzer::_PrintToFile\n"); 
#endif

	// lastly let's perform a little data output
	/***	This can be used as a header ***
	printf ("position");
	for (int atom = 0; atom < atomNames.size(); atom++) {
		printf ("%10s", atomNames[atom].c_str());
	} printf ("\n");
	*/

	// starting from the beginning of the file (i.e. overwrite it)
	rewind (_output);

	for (int i=0; i < _posbins; i++) {
		fprintf (_output, "% 10.4f", double(i)*_binsize+_start);		// the bin's position value

		// the output value of the number density will be converted to the actual density in g/mL.
		// the volume of a slice of the slab
		double volume = Atom::Size()[x] * Atom::Size()[y] * Atom::Size()[z] / Atom::Size()[_axis];
		volume *= _binsize;

		RUN2 (_atomNames) {
			// The density is thus transformed. The resulting values are the densities (mol/mL) of each species.
			// The value of 1.6611293 comes from a combo of avogadro's number and the conversion to mL from angstroms^3.
			// To get the density in g/mL, just multiply this value by the molecular/atomic weight of the species.
			// For molarity (mol/L) multiply this value by 1000 (to convert from mL to L)
			double density = double(_density[j][i]) / volume / ((double)_step+1.0) * 1.6611293;
			#ifdef AVG
			density /= 2.0;		// if we're averaging 2 interfaces, then these values need to be halved
			#endif

			fprintf (_output, "% 13.7f", density);			// for each atom printout the histogram value at that position
		}

		fprintf (_output, "\n");
	}

return;
}
	
// output a status meter to the screen
void DensityAnalyzer::_PrintStatus (int const step) {

#ifdef DEBUG
this->Debug ("DensityAnalyzer::_PrintStatus\n"); 
#endif

	if (!(step % (OUTPUT_FREQ * 10)))
		printf ("\n%10d/%d)  ", step, TIMESTEPS);

	if (!(step % OUTPUT_FREQ)) 
		printf ("*");

	fflush (stdout);

return;
}

vector<int> DensityAnalyzer::AtomDensity (string const atomname) {

#ifdef DEBUG
this->Debug ("DensityAnalyzer::AtomDensity\n"); 
#endif

	vector<int> density (_posbins, 0);
	
	// and now run through the actual calculations to find number densities
	RUN (_sys) {
			
		// find the atom of interest
		Atom * patom = _sys[i];
		
		//if (_sys[i]->Name().find(atomname) == string::npos) continue;
		if (_sys[i]->Name() != atomname) continue;

		// grab the position info of the atom
		VecR r = patom->Position();
		double position = r[_axis];
		if (position < PBCFLIP) position += Atom::Size()[_axis];

		#ifdef AVG
		// here the bin will be selected based on the distance to a given interface. Negative distances are inside the water phase, positive are in the CCl4
		if (position < START or position > END) continue;		// only bin stuff within the bounds that have been set up
		double distance = (position > middle) ? position - int_high : int_low - position;
		int bin = (int)((distance - _start)/_binsize);
		#else
		int bin = (int)((position - _start)/_binsize);
		#endif

		if (bin > _posbins) {
			printf ("DensityAnalyzer::AtomDensity - something is wrong - not enough bins in the histogram for analysis?\n");
			exit(1);
		}
		// and bin it into the density histogram
		density[bin]++;

	}

return (density);
}

void DensityAnalyzer::SystemDensities () {

#ifdef DEBUG
this->Debug ("DensityAnalyzer::SystemDensities\n"); 
#endif

	// now let's run through the timestesp
	for (_step=0; _step < _steps; _step++) {
		
		// for each atom that we're testing we'll add the histogram data into the final data-set
		for (unsigned int atom = 0; atom < _atomNames.size(); atom++) {
			vector<int> atomDensity = this->AtomDensity (_atomNames[atom]);

			// once we have the data for each atom for each timestep, let's add it into the running total
			for (unsigned int bin = 0; bin < _posbins; ++bin) {
				_density[atom][bin] = _density[atom][bin] + atomDensity[bin];
			}
		}
			
		// and set up the system for the next timestep
		_sys.LoadNext();
		this->_PrintStatus (_step);
		if (!(_step % (OUTPUT_FREQ * 25)))
			this->_PrintToFile ();

	}

	this->_PrintToFile ();
	fclose(_output);

return;
}


int main (int argc, char **argv) {

	DensityAnalyzer den (argv, argc, TIMESTEPS, START, END, BINSIZE, AXIS);

#ifdef DEBUG
den.Debug ("main\n");
#endif

	den.SystemDensities();

return 0;
}

void DensityAnalyzer::Debug (string msg) const {

	printf ("%s", msg.c_str());
	printf ("_density[0].size() == %d\n", _density[0].size());
	//printf ("&_density[0] == %d\n", &_density[0]);

return;
}
