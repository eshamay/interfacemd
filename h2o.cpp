#include "h2o.h"

int Water::numWaters = 0;

#ifdef H2O_DIPOLE_PARM
WaterDipoleParms Water::_dipparms ("dipoleparm.dat");
#endif

Water::Water ()
{
	this->Rename("h2o");
	_moltype = Molecule::H2O;
	++numWaters;
}

Water::~Water () {
	numWaters--;
}

Water::Water (const Molecule& mol) : Molecule(mol) {
	this->Rename("h2o");
	_moltype = Molecule::H2O;
	++numWaters;
}

Water::Water (const MolPtr& mol) : Molecule(*mol) {
	this->Rename("h2o");
	_moltype = Molecule::H2O;
	++numWaters;
}

void Water::SetBondLengths () {
	this->_oh1 = this->_h1->Position() - this->_o->Position();
	this->_oh2 = this->_h2->Position() - this->_o->Position();
	return;
}

void Water::SetAtoms () {

	if (!_set) {

		// first let's grab pointers to the three atoms and give them reasonable names
		this->_h1 = (AtomPtr)NULL; this->_h2 = (AtomPtr)NULL;

		for (Atom_it it = _atoms.begin(); it != _atoms.end(); it++) {
			if ((*it)->Element() == Atom::O)
				this->_o = *it;
			if ((*it)->Element() == Atom::H) {
				if (_h1 == (AtomPtr)NULL)
					this->_h1 = *it;
				else
					this->_h2 = *it;
			}

			if (*it == (AtomPtr)NULL) {
				std::cout << "problem setting the water atoms! Water::SetAtoms()" << std::endl;
			}
		}

		// we can calculate the two O-H vectors
		this->SetBondLengths ();

		this->Set();
	}
	return;
}

// flip the water about a given plane (perpendicular to the given axis) running through the oxygen
void Water::Flip (const coord axis) {

	this->SetAtoms();

	double center = _o->Position()[axis];
	double distance1 = _oh1[axis];
	double distance2 = _oh2[axis];

	VecR offset1, offset2;

	offset1.Set(axis, center - distance1);
	offset2.Set(axis, center - distance2);

	_h1->Position(axis, center - distance1);
	_h2->Position(axis, center - distance2);

	return;
}

VecR Water::Bisector () {

	this->SetAtoms();

	VecR bisector = _oh1.normalized() + _oh2.normalized();

	return bisector.normalized();
}

/* Another type of molecular axes (non-morita, useful for molecular symmetry/order parameter stuff) is as follows:
 * z-axis = the negative of the C2V (bisector) axis
 * y-axis = perpendicular to the plane of the molecule
 * x-axis = y % z
 */
void Water::SetOrderAxes () {

	this->SetAtoms ();

	// the z-axis is the negative of the C2V axis - so find the bisector and set the vector pointing towards the O (just like the dipole vector)
	_z = this->Bisector() * (-1.0);

	// the y-axis points perpendicular to the plane of the molecule. This can be found from the cross product of the two OH vectors
	_y = (_oh1 % _oh2).normalized();

	// and the x-axis is easy
	_x = (_y % _z).normalized();

	return;
}


VecR Water::MolecularAxis () {
	this->SetOrderAxes ();
	return _z;
}


