#ifndef ATOM_H_
#define ATOM_H_

#include <string>
#include <vector>
#include "vecr.h"
#include "utility.h"

class Molecule;

class Atom {

	string _name, 
		   _residue;
	
	int    _ID;				// some numerical identifier in case the atom is in an ordered list
	int	   _molid;		// the molecule that contains this atom

	Molecule * _pmolecule;

	double _mass, 
		   _charge;

	VecR _position;			// Particle position
	VecR _force;			// the external force on the atom at any given point in time

	static VecR _size;				// system size

public:
	// constructors
	Atom ();
	Atom (string name, VecR position);
	Atom (string name, VecR position, VecR force);
	Atom (VecR position);
	Atom (const Atom& oldAtom);				// copy constructor for deep copies

	double operator- (const Atom& input) const;		// operator usage to determine the distance between two atoms
	double operator[] (const coord index) const;	// get the atom's position by coordinate

	// Input
	void Name (const string name) { _name = name; }

	void Position (const VecR& position) { _position = position; }
	void Position (double X, double Y, double Z) { _position.Set(X, Y, Z); }
	void Position (coord const axis, double const value) { _position.Set (axis, value); }

	void Force (const VecR& force) { _force = force; }
	void Force (double X, double Y, double Z) { _force.Set(X, Y, Z); }
	void Force (coord const axis, double const value) { _force.Set (axis, value); }

	void ID (int id) { _ID = id; }
	//void Charge (double charge) { _charge = charge; }
	void SetCharge ();
	void Residue (string residue) { _residue = residue; }

	void X (double val) { _position.X(val); }			// for setting the atom's position
	void Y (double val) { _position.Y(val); }
	void Z (double val) { _position.Z(val); }
	void SetMass ();
	void MolID (const int mol) { _molid = mol; }	// sets the ID of the molecule containing this atom
	void ParentMolecule (Molecule * mol) { _pmolecule = mol; }	// sets a pointer to the molecule that contains the atom
	static void Size (VecR size) { _size = size; }	// sets the system size
	static VecR Size ()	{ return _size; }

	void Shift (VecR shift);			// shift the atom's position
	
	// Output
	string Name () const 	{ return (_name); }
	double Mass () const 	{ return _mass; }
	double Charge () const 	{ return _charge; }
	int ID () const 		{ return _ID; }
	string Residue () const { return _residue; }

	const VecR& Position () const	{ return _position; }
	//std::vector<double>& DPosition () { return _position.Coords(); }	// for returning the double array instead of the vector object

	const VecR& Force () const		{ return _force; }
	//double * DForce () 		{ return _force.Coords(); }

	double X () const 		{ return _position[x]; }
	double Y () const		{ return _position[y]; }
	double Z () const 		{ return _position[z]; }
	int MolID () const		{ return _molid; }
	Molecule * ParentMolecule () const { return _pmolecule; }
	void Wrap (VecR origin);							// A way to wrap the atom's position into the central image of the system cell
	void Print () const;
};

typedef std::vector<Atom>::iterator ATOM_IT;
typedef std::vector<Atom *>::iterator PATOM_IT;
typedef std::vector<Atom *> VPATOM;
typedef std::vector<Atom> VATOM;

#endif
