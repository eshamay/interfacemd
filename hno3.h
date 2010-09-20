#ifndef HNO3_H_
#define HNO3_H_

#include "vecr.h"
#include "molecule.h"
#include "mdsystem.h"
#include <iostream>


// A nitric acid class to add a few functions for dealing with nitric acid molecules specifically.
class NitricAcid: public Molecule {

	protected:
		// Nitric acid molecules have a molecular plane formed by the three oxygen molecules. This vector defines that plane
		VecR _molPlane;
		VecR _no2dipole;

		VecR_vec _no2wanniers;

		bool _set;								// to let us know if the atom have been set in the molecule
		AtomPtr _oh, _n, _h, _o1, _o2;		// Pointers to the various atoms in the nitric acid. These are static for no
		VecR _voh;
		// meaning that they're set at the the first formation of the molecule

	public:
		NitricAcid ();	// a default constructor
		~NitricAcid (); // deconstructor

		static int numNitricAcids;
		bool CalcNO2Dipole ();	// calculate the nitric acid dipole

		void SetAtoms ();
		AtomPtr GetOH () const { return _oh; }
		AtomPtr GetH () const { return _h; }

		// Functions for analysis
		VecR MolecularPlaneVector ();
		VecR MolecularAxis () { return _z; }
		VecR NO2Dipole () const { return _no2dipole; }
		VecR NO2Bisector ();
		void PrintNO2 () const;
		VecR& OH () { return _voh; }
		const VecR_vec& NO2Wanniers () const { return _no2wanniers; }

};


class Nitrate: public Molecule {
	protected:
		VecR _molPlane;

		bool _set;								// to let us know if the atom have been set in the molecule
		AtomPtr _n, _o1, _o2, _o3;		// Pointers to the various atoms in the nitric acid. These are static for no
		// meaning that they're set at the the first formation of the molecule

		VecR	_no1, _no2, _no3;		// vectors of the N-O bonds

	public:
		Nitrate ();	// a default constructor
		~Nitrate (); // deconstructor

		static int numNitrates;

		void SetAtoms ();

		VecR MolecularAxis () { return _z; }
		VecR const * NO1 () { return &_no1; }
		VecR const * NO2 () { return &_no2; }
		VecR const * NO3 () { return &_no3; }

		// Functions for analysis
		VecR MolecularPlaneVector ();

};


#endif
