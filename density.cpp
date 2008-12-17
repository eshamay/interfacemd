#include "density.h"

vector<int> NumberDensity (AmberSystem& sys, double const start, double const end, double const binsize, coord const axis) {

	// first we figure out how many bins there are on each axis
	int size = (end - start)/binsize + 1;

	// then create the output vector
	vector<int> density;

	// and we'll size and zero out the containers
	density.resize (size, 0);

	// and now run through the actual calculations to find number densities
	RUN (sys) {
			
		VecR r = sys[i]->Position();
		double position = r[axis];
		int bin = (int)((position-start)/binsize);

		density[bin]++;
	}

return (density);
}

vector<int> MoleculeDensity ( AmberSystem& sys, double const start, double const end, double const binsize, const coord axis, string atomname) {

	// first we figure out how many bins there are on each axis
	int size = (end - start)/binsize + 1;

	// then create the output vector
	vector<int> density;

	// and we'll size and zero out the containers
	density.resize (size, 0);

	// and now run through the actual calculations to find number densities
	RUN (sys.Molecules()) {
			
		Atom * patom = sys.Molecules(i)->GetAtom(atomname);
		VecR r = patom->Position();
		double position = r[axis];
		// temp stuff for a specific system
		if (position < 30.0) position += Atom::Size()[axis];
		int bin = (int)((position-start)/binsize);

		density[bin]++;
	}

return (density);
}
