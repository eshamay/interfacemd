#include "mdsystem.h"

// The system size for periodic boundary calculations
VecR MDSystem::_dimensions = VecR ();

MDSystem::~MDSystem () {
  return;
}

// Find the smallest vector between two locations in a periodic system defined by the dimensions.
// The resulting vector will point from the v1 to v2
VecR MDSystem::Distance (const VecR& v1, const VecR& v2) {

	// first we gather all our coordinates for point a (the starting point) and point b (the ending point)
	double ax = v1[x];
	double ay = v1[y];
	double az = v1[z];

	double bx = v2[x];
	double by = v2[y];
	double bz = v2[z];

	// Now we'll hold one coordinate while moving the other through its periodic images until the distance between the two is <= size/2
	// this is very much like the fmod() function, but it's written out for clarity here to show that one is being held fixed
	while (fabs(ax-bx) > _dimensions[x]/2.0) {
		if (ax < bx) ax += _dimensions[x];
		else 		 ax -= _dimensions[x];
	}

	while (fabs(ay-by) > _dimensions[y]/2.0) {
		if (ay < by) ay += _dimensions[y];
		else 		 ay -= _dimensions[y];
	}

	while (fabs(az-bz) > _dimensions[z]/2.0) {
		if (az < bz) az += _dimensions[z];
		else 		 az -= _dimensions[z];
	}

return (VecR(bx - ax, by - ay, bz - az));
}

// A useful method for finding the distances between pairs of atoms in a periodic system
VecR MDSystem::Distance (const Atom * atom1, const Atom * atom2) {
	return MDSystem::Distance(atom1->Position(), atom2->Position());
}


