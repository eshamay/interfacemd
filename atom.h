#pragma once
#ifndef ATOM_H_
#define ATOM_H_

#include "vecr.h"
#include "utility.h"
#include <vector>
#include <string>

class Molecule;
typedef Molecule* MolPtr;


class Atom {

  public:

	typedef enum {
	  NO_ELEMENT = 0,
	  H = 1, He = 2,
	  B = 5, C = 6, N = 7, O = 8, F = 9, Ne = 10,
	  Na = 11, Mg = 12, Al = 13, Si = 14, P = 15, S = 16, Cl = 17, Ar = 18,
	  K = 19, Ca = 20, 
	  I = 53, 
	  Cs = 55, Ba = 56,
	  D, DW, SW
	} Element_t;

	// constructors
	Atom ();
	Atom (const std::string& name, const VecR& position);
	Atom (const std::string& name, const VecR& position, const VecR& force);
	Atom (const VecR& position);
	Atom (const Atom& oldAtom);				// copy constructor for deep copies
	virtual ~Atom ();

	typedef Atom* AtomPtr;
	typedef std::vector<AtomPtr> Atom_ptr_vec;
	typedef Atom_ptr_vec::const_iterator Atom_it;

	double operator- (const Atom& input) const;		// operator usage to determine the distance between two atoms
	double operator[] (const coord index) const;	// get the atom's position by coordinate
	bool operator< (const AtomPtr& rhs) const;
	bool operator< (const Atom& rhs) const;

	// Input
	virtual void UnsetMolecule() const;

	void Name (const std::string& name) { _name = name; UnsetMolecule(); }

	void Position (const VecR& position);
	void Position (const double X, const double Y, const double Z);

	void Position (coord const axis, double const value);

	void Force (const VecR& force) { _force = force; UnsetMolecule(); }
	void Force (const double X, const double Y, const double Z) { _force.Set(X, Y, Z); UnsetMolecule(); }
	void Force (coord const axis, double const value) { _force.Set (axis, value); UnsetMolecule(); }

	void ID (const int id) { _ID = id; UnsetMolecule(); }
	//void Charge (double charge) { _charge = charge; }
	void SetAtomProperties ();
	void Residue (const std::string& residue) { _residue = residue; UnsetMolecule(); }

	// for setting the atom's position
	void X (double val) { _position[0] = val; UnsetMolecule(); }			
	void Y (double val) { _position[1] = val; UnsetMolecule(); }
	void Z (double val) { _position[2] = val; UnsetMolecule(); }

	void MolID (const int mol) { _molid = mol; UnsetMolecule(); }	// sets the ID of the molecule containing this atom
	void ParentMolecule (const MolPtr mol) { _pmolecule = mol; UnsetMolecule(); }	// sets a pointer to the molecule that contains the atom

	void Shift (const VecR& shift)			// shift the atom's position
	{ _position += shift; UnsetMolecule(); }

	// Output
	const std::string& Name () const 	{ return _name; }
	const Element_t& Element () const { return _element; }
	const double& Mass () const 	{ return _mass; }
	const double& Charge () const 	{ return _charge; }
	const int& ID () const 		{ return _ID; }
	const std::string& Residue () const { return _residue; }

	const VecR& Position () const	{ return _position; }

	const VecR& Force () const		{ return _force; }

	double X () const 		{ return _position.x(); }
	double Y () const		{ return _position.y(); }
	double Z () const 		{ return _position.z(); }
	int MolID () const		{ return _molid; }
	MolPtr ParentMolecule () const { return _pmolecule; }
	void Print () const;


	static AtomPtr FindByElement (const Atom_ptr_vec& apv, Element_t elmt) {
		Atom_it a = std::find_if (apv.begin(), apv.end(), member_functional::mem_fun_eq(&Atom::Element,elmt));
		return *a;
	}

	static bool element_eq (const AtomPtr& first, const AtomPtr& second) {
		return first->Element() == second->Element();
	}

	// tests if the combination of atoms supplied matches the element pair supplied

	// some predicates
	static bool ElementCombo (const AtomPtr& ai, const AtomPtr& aj, const Element_t element_a, const Element_t element_b) {
		return 
			((ai->Element() == element_a && aj->Element() == element_b)
			 ||
			 (ai->Element() == element_b && aj->Element() == element_a));
	}



	protected:
	std::string _name, 	// human-readable identifier
		_residue; 		// name of the parent-molecule 


	int    _ID;	// some numerical identifier in case the atom is in an ordered list
	int	   _molid;			   // the molecule that contains this atom

	MolPtr _pmolecule;

	double _mass, _charge;

	Element_t _element;			// the actual element based on the atom name - always upper-case and max length of two letters

	VecR _position;				// Particle position
	VecR _force; // the external force on the atom at any given point in time

};	// class Atom


typedef Atom::AtomPtr AtomPtr;
typedef Atom::Atom_ptr_vec Atom_ptr_vec;
typedef Atom::Atom_it Atom_it;

#endif
