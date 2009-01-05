#include "densitytest.h"

// the command line arguments are the list of atom *NAMES* that are to be analyzed... exactly as they are named in the prmtop file
// the files prmtop, mdcrd, and mdvel should point to their respective files.
DensityAnalyzer::DensityAnalyzer (char * argv[], int const argc, int const numSteps, double const start, double const end, double const binsize, coord const axis) : _sys(AmberSystem(PRMTOP, MDCRD, FORCE)) {

	// this will be our output data file
	_output = (FILE *)NULL;
	_output = fopen ("density.dat", "w");
	if (_output == (FILE *)NULL) {
		printf ("couldn't open the output file density.dat. Now Exiting\n");
		exit (1);
	}

	_step		= 0;
	_steps 		= numSteps;
	_start 		= start;
	_end 		= end;
	_binsize 	= binsize;
	_axis		= axis;

	printf ("Performing a DENSITY analysis of the system\n");
	printf ("Position bounds are set at: % 8.3f to % 8.3f in % 8.3f increments\n", _start, _end, _binsize);

	// first we figure out how many bins there are on each axis
	_posbins = int ((_end - _start)/_binsize) + 1;

	// the location of the two interfaces
	int_low = INT_LOW;
	int_high = INT_HIGH;
	middle = (int_low + int_high)/2.0;

	printf ("Interfaces are set at % 8.3f and % 8.3f\n", int_low, int_high);
	
	// from the command line, grab all the atom name that we'll be working with
	for (int i = 1; i < argc; i++) {
		_atomNames.push_back (argv[i]);
	}

	// now we set up the histogram(s) such that we may bin the positions of each atom
	_density.resize(_atomNames.size(), vector<int> (_posbins, 0));
	
return;
}

void DensityAnalyzer::_PrintToFile () {

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
			// The density is thus transformed. 
			// The value of 1.6611293 comes from a combo of avogadro's number and the conversion to mL from angstroms^3.
			// To get the density in g/mL, just divide this value by the molecular/atomic weight of the species.
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

	if (!(step % (OUTPUT_FREQ * 10)))
		printf ("\n%10d/%d)  ", step, TIMESTEPS);

	if (!(step % OUTPUT_FREQ)) 
		printf ("*");
	
	fflush (stdout);

return;
}

vector<int> DensityAnalyzer::AtomDensity (string const atomname) {

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
		if (position < 17.0) position += Atom::Size()[_axis];

		// here the bin will be selected based on the distance to a given interface. Negative distances are inside the water phase, positive are in the CCl4
		double distance = (position > middle) ? position - int_high : int_low - position;
		
		// and bin it into the density histogram
		#ifdef AVG
		int bin = (int)((distance - _start)/_binsize);
		#else
		int bin = (int)((position - _start)/_binsize);
		#endif

		density[bin]++;
	}
return (density);
}

void DensityAnalyzer::SystemDensities () {

	// some info before starting
	printf ("Analyzing densities of the atoms: ");
	RUN (_atomNames) {
		printf ("%s  ", _atomNames[i].c_str());
	}
	printf ("\nOutputting data to file: %s", "density.dat");


	// now let's run through the timestesp
	for (_step=0; _step < _steps; _step++) {
		
		// for each atom that we're testing we'll add the histogram data into the final data-set
		RUN2 (_atomNames) {
			vector<int> atomDensity = this->AtomDensity (_atomNames[j]);

			// once we have the data for each atom for each timestep, let's add it into the running total
			for (int i=0; i < _posbins; i++) {
				_density[j][i] += atomDensity[i];
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

	den.SystemDensities();

return 0;
}
