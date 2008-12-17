#include "densitytest.h"

// the command line arguments are the list of atom *NAMES* that are to be analyzed... exactly as they are named in the prmtop file
// the files prmtop, mdcrd, and mdvel should point to their respective files.
DensityAnalyzer::DensityAnalyzer (char * argv[], int const argc, int const numSteps, double const start, double const end, double const binsize, coord const axis) : _sys(AmberSystem(PRMTOP, MDCRD, FORCE)) {

	// this will be our output data file
	_output = fopen ("density.dat", "w");

	_steps 		= numSteps;
	_start 		= start;
	_end 		= end;
	_binsize 	= binsize;
	_axis		= axis;

	// first we figure out how many bins there are on each axis
	_size = int ((_end - _start)/_binsize) + 1;
	
	// from the command line, grab all the atom name that we'll be working with
	for (int i = 1; i < argc; i++) {
		_atomNames.push_back (argv[i]);
	}

	// now we set up the histogram(s) such that we may bin the positions of each atom
	_density.resize(_atomNames.size(), vector<int> (_size, 0));
	
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

	for (int i=0; i < _size; i++) {
		fprintf (_output, "% 10.4f", double(i)*_binsize+_start);		// the bin's position value

		for (unsigned int atom = 0; atom < _atomNames.size(); atom++) {
			fprintf (_output, "% 10d", _density[atom][i]);			// for each atom printout the histogram value at that position
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

	vector<int> density (_size, 0);
	
	// and now run through the actual calculations to find number densities
	RUN (_sys) {
			
		// find the atom of interest
		Atom * patom = _sys[i];
		//if (_sys[i]->Name().find(atomname) == string::npos) continue;
		if (_sys[i]->Name() != atomname) continue;

		// grab the position info of the atom
		VecR r = patom->Position();
		double position = r[_axis];
		if (position < 15.0) position += Atom::Size()[_axis];

		// and bin it into the density histogram
		int bin = (int)((position - _start)/_binsize);

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
	for (int step=0; step < _steps; step++) {
		
		// for each atom that we're testing we'll add the histogram data into the final data-set
		for (unsigned int atom=0; atom < _atomNames.size(); atom++) {
			vector<int> atomDensity = this->AtomDensity (_atomNames[atom]);

			// once we have the data for each atom for each timestep, let's add it into the running total
			for (int i=0; i < _size; i++) 
				_density[atom][i] += atomDensity[i];
		}
			
		// and set up the system for the next timestep
		_sys.LoadNext();

		this->_PrintStatus (step);
		if (!(step % (OUTPUT_FREQ * 25)))
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
