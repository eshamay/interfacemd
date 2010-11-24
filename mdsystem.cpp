#include "mdsystem.h"

// The system size for periodic boundary calculations
VecR MDSystem::_dimensions = VecR ();

MDSystem::~MDSystem () {
	return;
}

// Find the smallest vector between two locations in a periodic system defined by the dimensions.
// The resulting vector will point from the v1 to v2
VecR MDSystem::Distance (const VecR v1, const VecR v2) {

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
VecR MDSystem::Distance (const AtomPtr atom1, const AtomPtr atom2) {
	return MDSystem::Distance(atom1->Position(), atom2->Position());
}

double MDSystem::Distance (const MolPtr mol1, const MolPtr mol2) {

	std::vector <double> distances;
	for (Atom_it ai = mol1->begin(); ai != mol1->end(); ai++) {
		for (Atom_it aj = mol2->begin(); aj != mol2->end(); aj++) {
			distances.push_back (MDSystem::Distance(*ai,*aj).Magnitude());
		}
	}
	std::sort(distances.begin(), distances.end());
	return distances[0];
}

VecR MDSystem::CalcClassicDipole (MolPtr mol) {
	VecR com = mol->UpdateCenterOfMass();
	VecR dipole = VecR::Zero();

	// the dipole is just a sum of the position vectors multiplied by the charges (classical treatment)
	for (Atom_it it = mol->begin(); it != mol->end(); it++) {
		VecR r (MDSystem::Distance(com, (*it)->Position()));
		r *= (*it)->Charge();
		//printf ("%4.3f  % 4.3f ", r.Magnitude(), (*it)->Charge()); r.Print();
		dipole += r;
	}

	mol->Dipole(dipole);

	return (dipole);
}

VecR MDSystem::CalcWannierDipole (MolPtr mol) {

	VecR dipole = CalcClassicDipole(mol);

	VecR com = mol->CenterOfMass();

	// wannier centers have a charge of -2
	for (VecR_it vt = mol->wanniers_begin(); vt != mol->wanniers_end(); vt++) {
		VecR r (MDSystem::Distance(com, *vt));
		r *= 2.0;
		dipole -= r;
	}

	//mol->Print();
	//printf ("% 8.3f ) ", dipole.Magnitude()); dipole.Print();
	mol->Dipole(dipole);

	return (dipole);
}


