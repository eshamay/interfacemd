#ifndef ATOM_H_
#define ATOM_H_

#include <string>
#include <vector>
#include "vecr.h"

class Molecule;
typedef Molecule* MolPtr;

class Atom {

  protected:
    std::string _name,
      _residue;

    int    _ID;	// some numerical identifier in case the atom is in an ordered list
    int	   _molid;			   // the molecule that contains this atom

    MolPtr _pmolecule;

    double _mass,
	   _charge;

    VecR _position;				// Particle position
    VecR _force; // the external force on the atom at any given point in time

  public:

    // constructors
    Atom ();
    Atom (std::string name, VecR position);
    Atom (std::string name, VecR position, VecR force);
    Atom (VecR position);
    Atom (const Atom& oldAtom);				// copy constructor for deep copies
    ~Atom ();

	typedef Atom * AtomPtr;

    double operator- (const Atom& input) const;		// operator usage to determine the distance between two atoms
    double operator[] (const coord index) const;	// get the atom's position by coordinate

    // Input
    void Name (const std::string name) { _name = name; }

    void Position (const VecR& position) { _position = position; }
    void Position (double X, double Y, double Z) { _position.Set(X, Y, Z); }

    void Position (coord const axis, double const value) { _position.Set (axis, value); }

    void Force (const VecR& force) { _force = force; }
    void Force (double X, double Y, double Z) { _force.Set(X, Y, Z); }
    void Force (coord const axis, double const value) { _force.Set (axis, value); }

    void ID (int id) { _ID = id; }
    //void Charge (double charge) { _charge = charge; }
    void SetCharge ();
    void Residue (std::string residue) { _residue = residue; }

    // for setting the atom's position
    void X (double val) { _position.X(val); }			
    void Y (double val) { _position.Y(val); }
    void Z (double val) { _position.Z(val); }

    void SetMass ();
    void MolID (const int mol) { _molid = mol; }	// sets the ID of the molecule containing this atom
    void ParentMolecule (Molecule * mol) { _pmolecule = mol; }	// sets a pointer to the molecule that contains the atom

    void Shift (VecR shift)			// shift the atom's position
    {
      _position += shift;
    }

    // Output
    std::string Name () const 	{ return (_name); }
    double Mass () const 	{ return _mass; }
    double Charge () const 	{ return _charge; }
    int ID () const 		{ return _ID; }
    std::string Residue () const { return _residue; }

    const VecR& Position () const	{ return _position; }

    const VecR& Force () const		{ return _force; }

    double X () const 		{ return _position[x]; }
    double Y () const		{ return _position[y]; }
    double Z () const 		{ return _position[z]; }
    int MolID () const		{ return _molid; }
    MolPtr ParentMolecule () const { return _pmolecule; }
    void Print () const;

	class NameIs_p : public std::binary_function<AtomPtr,std::string,bool> {
	  public:
		bool operator() (const AtomPtr atom, const std::string name) const {
		  return atom->Name() == name;
		}
	};

	class AtomPtr_sort : public std::binary_function<AtomPtr,AtomPtr,bool> {
	  public:
		bool operator() (AtomPtr const &left, AtomPtr const &right) { 
		  return left->ID() < right->ID();
		}
	};

};

typedef Atom::AtomPtr AtomPtr;
typedef std::vector<AtomPtr> Atom_ptr_vec;
typedef Atom_ptr_vec::const_iterator Atom_it;

#endif
