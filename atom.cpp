#include "atom.h"

VecR Atom::_size (0.0, 0.0, 0.0);

Atom::Atom () :
	_name(""),
	_residue(""),
	_ID(-1),
	_molid(-1),
	_pmolecule((Molecule *)NULL),
	_mass(0.0),
	_charge(0.0),
	_position(VecR()),
	_force(VecR()) {
}

Atom::Atom (std::string name, VecR position, VecR force) :
	_name(name),
	_residue(""),
	_ID(-1),
	_molid(-1),
	_pmolecule((Molecule *)NULL),
	_mass(0.0),
	_charge(0.0),
	_position(position),
	_force(force) {

	this->SetMass ();			// set the atom's mass if the name is known
	this->SetCharge ();
}

Atom::Atom (std::string name, VecR position) :
	_name(name),
	_residue(""),
	_ID(-1),
	_molid(-1),
	_pmolecule((Molecule *)NULL),
	_mass(0.0),
	_charge(0.0),
	_position(position),
	_force(VecR()) {

	this->SetMass ();			// set the atom's mass if the name is known
	this->SetCharge ();
}

// copy constructor will pull over all the member values, and copy over the vector values
Atom::Atom (const Atom& oldAtom) :
	_name(oldAtom._name),
	_residue(oldAtom._residue),
	_ID(oldAtom._ID),
	_molid(oldAtom._molid),
	_pmolecule(oldAtom._pmolecule),
	_mass(oldAtom._mass),
	_charge(oldAtom._charge),
	_position(oldAtom._position),
	_force(oldAtom._force) {

	this->SetMass();
	this->SetCharge();
}

Atom::Atom (VecR position) :
	_name(""),
	_residue(""),
	_ID(-1),
	_molid(-1),
	_pmolecule((Molecule *)NULL),
	_mass(0.0),
	_charge(0.0),
	_position(position),
	_force(VecR()) {
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

double Atom::MinDistance (const Atom& input) const {
	return (_position.MinDistance(input.Position(), _size));
}

double Atom::MinDistance (Atom const * const input) const {
	return (_position.MinDistance(input->Position(), _size));
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

