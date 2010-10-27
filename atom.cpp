#include "atom.h"

Atom::Atom () :
  _name(""), _residue(""), _ID(-1), _molid(-1),
  _pmolecule((Molecule *)NULL),
  _mass(0.0), _charge(0.0), _element(NO_ELEMENT)
{ }

Atom::Atom (const std::string& name, const VecR& position, const VecR& force) :
  _name(name), _residue(""), _ID(-1), _molid(-1),
  _pmolecule((Molecule *)NULL),
  _mass(0.0), _charge(0.0),
  _position(position), _force(force) 
{ this->SetAtomProperties (); }

Atom::Atom (const std::string& name, const VecR& position) :
  _name(name),
  _residue(""),
  _ID(-1),
  _molid(-1),
  _pmolecule((Molecule *)NULL),
  _mass(0.0),
  _charge(0.0),
  _position(position)
{ this->SetAtomProperties (); }

// copy constructor will pull over all the member values, and copy over the vector values
Atom::Atom (const Atom& oldAtom) :
  _name(oldAtom._name), _residue(oldAtom._residue), _ID(oldAtom._ID), _molid(oldAtom._molid),
  _pmolecule(oldAtom._pmolecule),
  _mass(oldAtom._mass), _charge(oldAtom._charge),
  _position(oldAtom._position), _force(oldAtom._force) 
{ this->SetAtomProperties (); }


Atom::Atom (const VecR& position) :
  _name(""), _residue(""), _ID(-1), _molid(-1),
  _pmolecule((Molecule *)NULL),
  _mass(0.0), _charge(0.0),
  _position(position)
{ }


Atom::~Atom () {
}

void Atom::Position (const VecR& position) { 
	_position = position; 
	UnsetMolecule();
}

void Atom::Position (const double X, const double Y, const double Z) { 
	_position.Set(X, Y, Z); 
	UnsetMolecule();
}

void Atom::Position (coord const axis, double const value) { 
	_position.Set (axis, value); 
	UnsetMolecule();
}

double Atom::operator[] (const coord index) const {
	double pos = -1.0;

	switch (index) {
		case x: pos = _position.x();
		case y: pos = _position.y();
		case z: pos = _position.z();
	}
	UnsetMolecule();
	return pos;
}

double Atom::operator- (const Atom& input) const {
	return (_position - input.Position()).Magnitude();
}

bool Atom::operator< (const AtomPtr& rhs) const {
	return this->_ID < rhs->ID();
}
bool Atom::operator< (const Atom& rhs) const {
	return this->_ID < rhs.ID();
}

void Atom::Print () const {
	printf ("%5s (ID:%6d nuc:%4d)\t%4s (molID:%5d)\t% f\t% f\t% f\n", _name.c_str(), _ID, _element, _residue.c_str(), _molid, _position.x(), _position.y(), _position.z());
}

void Atom::SetAtomProperties () {

	//const double AMU2KG	= 1.6762158e-27;		// conversion for amu -> kg
	if (_name.find("DW") != std::string::npos) {
		_mass = 0.0;
		_charge = 0.0;
		_element = DW;
	}
	else if (_name.find("SW") != std::string::npos) {
		_mass = 0.0;
		_charge = 0.0;
		_element = SW;
	}
	else if (_name.find("Cl") != std::string::npos) {
		_mass = 35.453;
		_charge = 7.0;
		_element = Cl;
	}
	else if (_name.find("H") != std::string::npos) {
		_mass = 1.0078250321;
		_charge = 1.0;
		_element = H;
	}
	else if (_name.find("O") != std::string::npos) {
		_mass = 15.9949146221;
		_charge = 6.0;
		_element = O;
	}
	else if (_name.find("N") != std::string::npos) {
		_mass = 14.0030740052;
		_charge = 5.0;
		_element = N;
	}
	else if (_name.find("D") != std::string::npos) {
		_mass = 2.0156500641;
		_charge = 1.0;
		_element = D;
	}
	else if (_name.find("C") != std::string::npos) {
		_mass = 12.0000000;
		_charge = 4.0;
		_element = C;
	}
	else if (_name.find("Si") != std::string::npos) {
		_mass = 28.0855;
		_charge = 4.0;
		_element = Si;
	}
	else if (_name.find("S") != std::string::npos) {
		_mass = 32.065;
		_charge = 6.0;
		_element = S;
	}
	else if (_name.find("F") != std::string::npos) {
		_mass = 18.9984;
		_charge = 7.0;
		_element = F;
	}

	UnsetMolecule();

	return;
}

