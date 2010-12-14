#include "so2.h"

SulfurDioxide::SulfurDioxide () {
	this->Rename("so2");
	_moltype = Molecule::SO2;
}


SulfurDioxide::SulfurDioxide (const Molecule& mol)
	:
		Molecule (mol)
{
	this->Rename("so2");
	_moltype = Molecule::SO2;
}


SulfurDioxide::SulfurDioxide (const MolPtr& mol)
	:
		Molecule (*mol)
{
	this->Rename("so2");
	_moltype = Molecule::SO2;
}

void SulfurDioxide::SetAtoms () {

	_s = this->GetAtom(Atom::S);
	_o1 = (AtomPtr)NULL; _o2 = (AtomPtr)NULL;

	for (Atom_it it = this->begin(); it != this->end(); it++) {
		if ((*it)->Element() == Atom::O) {
			if (_o1 == NULL)
				_o1 = *it;
			else
				_o2 = *it;
		}
	}

	this->_so1 = this->_o1->Position() - this->_s->Position();
	this->_so2 = this->_o2->Position() - this->_s->Position();
	return;
}

VecR SulfurDioxide::Bisector () {

	this->SetAtoms();

	VecR bisector = _so1.normalized() + _so2.normalized();

	return bisector.normalized();
}


/* Setting molecular axes (useful for molecular symmetry/order parameter stuff) is as follows:
 * z-axis = the C2V (bisector) axis - points from the S to the Os
 * y-axis = perpendicular to the plane of the molecule
 * x-axis = y % z
 */
void SulfurDioxide::SetBisectorAxes () {

	this->SetAtoms ();

	// the z-axis is the negative of the C2V axis - so find the bisector and set the vector pointing towards the O (just like the dipole vector)
	_z = this->Bisector();

	// the y-axis points perpendicular to the plane of the molecule. This can be found from the cross product of the two OH vectors
	_y = (_so1 % _so2).normalized();

	// and the x-axis is easy
	_x = (_y % _z).normalized();

	return;
}

