#include "atom.h"

VecR Atom::_size (0.0, 0.0, 0.0);

Atom::Atom () {

	_name = "";
	_residue = "";
	_ID = -1;
	_molid = -1;
	_pmolecule = (Molecule *)NULL;
	_mass = 0.0;
	_charge = 0.0;

	_position = VecR();
	_force = VecR();
}

Atom::Atom (std::string name, VecR position, VecR force) {
	_name = name;
	_position = position;

	_residue = "";
	_ID = -1;
	_molid = -1;
	_pmolecule = (Molecule *)NULL;
	_charge = 0.0;
	_mass = 0.0;

	_force = force;

	this->SetMass ();			// set the atom's mass if the name is known
	this->SetCharge ();
}

Atom::Atom (std::string name, VecR position) {
	_name = name;
	_position = position;

	_residue = "";
	_ID = -1;
	_molid = -1;
	_pmolecule = (Molecule *)NULL;
	_charge = 0.0;
	_mass = 0.0;

	_force = VecR();

	this->SetMass ();			// set the atom's mass if the name is known
	this->SetCharge ();
}

// copy constructor will pull over all the member values, and copy over the vector values
Atom::Atom (const Atom& oldAtom) {
	_name = oldAtom.Name();
	_residue = oldAtom.Residue();
	_ID = oldAtom.ID();
	this->SetMass();
	this->SetCharge();

	// copying the vectors over
	_position = oldAtom.Position();
	_force = oldAtom.Force();
}

Atom::Atom (VecR position) {
	_position = position;
}

double Atom::operator[] (const coord index) const {
	double pos;

	switch (index) {
		case x: pos = _position[x];
		case y: pos = _position[y];
		case z: pos = _position[z];
	}
	return pos;
}

double Atom::operator- (const Atom& input) const {
	return (_position.MinDistance(input.Position(), _size));
}

void Atom::Print () const {
	printf ("%s (%d)\t%s\t% f\t% f\t% f\n", _name.c_str(), _ID, _residue.c_str(), _position[x], _position[y], _position[z]);
}

void Atom::SetMass () {

	//const double AMU2KG	= 1.6762158e-27;		// conversion for amu -> kg
	if (_name.find("O") != std::string::npos) _mass = 15.9949146221;
	if (_name.find("N") != std::string::npos) _mass = 14.0030740052;
	if (_name.find("H") != std::string::npos) _mass = 1.0078250321;
	if (_name.find("D") != std::string::npos) _mass = 2.0156500641;
	if (_name.find("C") != std::string::npos) _mass = 12.0000000;
	if (_name.find("S") != std::string::npos) _mass = 32.065;

return;
}

// wrap the atom's position to the central periodic-image, given by the size of the system
void Atom::Wrap (VecR origin = VecR(0.0, 0.0, 0.0)) {
	_position.Wrap(_size, origin);
return;
}

void Atom::Shift (VecR shift) {
	_position += shift;
}

void Atom::SetCharge () {
	
	if (_name.find("O") != std::string::npos) _charge =  6.0;
	if (_name.find("H") != std::string::npos) _charge = 1.0;
	if (_name.find("N") != std::string::npos) _charge =  5.0;
	if (_name.find("S") != std::string::npos) _charge =  6.0;
	if (_name.find("Cl") != std::string::npos) _charge =  7.0;
	if (_name.find("C") != std::string::npos) _charge =  4.0;

return;
}

