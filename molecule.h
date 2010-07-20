#ifndef MOLECULE_H_
#define MOLECULE_H_

#include "vecr.h"
#include "matrixr.h"
#include "atom.h"
#include "wannier.h"
#include <string>
#include <math.h>
#include <vector>

class Molecule {

  protected:
    Atom_ptr_vec	_atoms;				// the list of the atoms in the molecule
    VecR_vec		_wanniers;			// the wannier centers in the molecule
    VecR			_dipole;			// the molecular dipole
    VecR			_x, _y, _z;			// molecular frame axes

    bool			_set;				// just a little helper to see if the atoms of the molecule have been set or for any other special purpose

    // this is broken last I checked - not updated with coordinate updates
    VecR			_centerofmass;		// calculate by 1/M * Sum(m[i]*r[i])	where M = total mass, m[i] and r[i] are atom mass and pos

    double			_mass;				// Total molecular mass
    double			_charge;
	std::string			_name;				// some text ID or name for the molecule
    int				_ID;				// A numerical ID for the molecule

    double			_eulerangles[3];	// the three euler angles theta, phi, chi
    MatR			_DCM;				// the direction cosine matrix for rotating the molecule to the lab frame
    void 			_FindEulerAngles ();// Calculates the Euler angles between the molecular axes and the fixed axes


  public:

    Molecule ();								// an empty molecule
    Molecule (const Molecule& oldMol);			// A copy constructor for performing deep copies of a molecule
    virtual ~Molecule ();

	typedef Molecule* MolPtr;

    static int numMolecules;

    // Input functions
    void Name (std::string name) { _name = name; }	// set the molecule's name
    void MolID (int ID) { _ID = ID; }

    bool Set () { _set = true; return (_set); }
    bool Unset () { _set = false; return (_set); }
	virtual void SetAtoms ();

    // Controls
    void Shift (VecR& shift);				// Shift the origin of the entire molecule
    void clear ();							// Erases the molecule data
    VecR UpdateCenterOfMass ();				// recalculates the center of mass when coordinates are updated

    // Output Functions
    VecR CenterOfMass () const		{ return _centerofmass; }
    VecR Position () const			{ return _centerofmass; }

    /* Dealing with atoms in the molecule */
    Atom_ptr_vec Atoms () const			{ return _atoms; }
    //AtomPtr Atoms (int index) const		{ return _atoms[index]; }

    Atom_it begin() const {
      return _atoms.begin();
    }
    Atom_it end() const {
      return _atoms.end();
    }

    const VecR_vec& Wanniers ()		const { return _wanniers; }
	VecR_it wanniers_begin () const { return _wanniers.begin(); }
	VecR_it wanniers_end () const { return _wanniers.end(); }

    double Mass () const 			{ return _mass; }					// Returns the molecular mass
    int size () const				{ return _atoms.size(); }
	std::string Name () const			{ return _name; }
    int MolID () const				{ return _ID; }

    // molecular axes
    VecR X () const					{ return _x; }
    VecR Y () const					{ return _y; }
    VecR Z () const					{ return _z; }
    // setting molecular axes
    void X (VecR& x_axis) { _x = x_axis; }
    void Y (VecR& y_axis) { _y = y_axis; }
    void Z (VecR& z_axis) { _z = z_axis; }

    // Euler angles
    double * EulerAngles () 		{ return _eulerangles; }

    // Initially written to output all the atom position data as a single contiguous double array for MPI xfer
    //double * DPositions () const;
    //double * DForces () const;

    void Print () const;							// print out all the info of the molecule

    double MinDistance (Molecule& mol);	// calculates the minimum distance between this molecule and another (2 closest atoms)

    //VecR CalcDipole ();	// calculate the dipole
	virtual void Dipole (VecR& dip) { _dipole = dip; }
    virtual VecR Dipole () const { return _dipole; }		// return the dipole of the molecule

    virtual VecR MolecularAxis () { return _z; }

    // Operators
    AtomPtr operator[] (const int index) const { return _atoms[index]; }	// retrieve an atom by array index
    AtomPtr operator[] (const std::string& atomname) const;			// retrieve a particular atom using its unique name/ID
	AtomPtr operator[] (const Atom::Element_t elmt) const;
    AtomPtr GetAtom (const std::string& atomname) const;
	AtomPtr GetAtom (const Atom::Element_t elmt) const;
    //int operator+= (Atom * newAtom);					// adds an atom into the molecule

    void AddAtom (AtomPtr const newAtom);					// same as the operator
    void AddHydrogen (AtomPtr const atom);					// same as adding an atom but renames accordingly
    void RemoveAtom (AtomPtr const atom);
    void FixAtom (AtomPtr const atom);
	void Rename (const std::string& name);

	MolPtr Merge (MolPtr mol);				// merges two molecules
	//int operator+= (Molecule& mol);					// Joins two molecules

	// Some stuff to work with wannier centers
	void AddWannier (const VecR& wannier) { _wanniers.push_back(wannier); } // adds a wannier center into the molecule
	void ClearWanniers () { _wanniers.clear(); }	// clear out the entire list

	//void ClearHBonds ();
	// return all the Hbonds that this molecule is involved in
	//std::vector<Atom *> HBonds () const;
	//int NumHBonds () const { return this->HBonds().size(); }


	// some functions to manipulate the molecule's position/orientation (symmetry operations)
	void Reflect (coord const axis, double const plane = 0.0);
	void Rotate (VecR& origin, VecR& axis, double angle);

	// A couple cool things used to manipulate matrices and vectors (don't know why this is in here, but hell, it works for now)
	// Constructs a rotation matrix (direction cosince matrix) from given euler angles to rotate *FROM THE LAB FRAME* to the molecule frame
	void RotateToLab (double vector[3]) const;
	void RotateToLab (double vector[][3]) const;

	// Constructs a rotation matrix to rotate *TO THE LAB FRAME* from the molecule-fixed frame
	void RotateToMol (double vector[3]) const;
	void RotateToMol (double matrix[][3]) const;

	// Performs matrix (3x3) multiplication (rotation) on a vector (3x1)
	void RotateVector (double rotation[][3], double vector[]) const;

	// Performs matrix rotation (i.e. matrix multiplication)
	void RotateMatrix (double rotation[][3], double matrix[][3]) const;

	// get the rotation matrix to rotate a molecule to lab-frame coordinates
	MatR const & DCMToLab (const coord axis = z);
};

typedef Molecule::MolPtr MolPtr;
typedef std::vector<MolPtr> Mol_ptr_vec;
typedef Mol_ptr_vec::const_iterator Mol_it;

#endif
