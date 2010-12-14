#ifndef CARBONCHAIN_H_
#define CARBONCHAIN_H_

#include "molecule.h"

	class Alkane : public Molecule {

		protected:
			Atom_ptr_vec _carbons;			/* An ordered listing of all the carbons in
																		 the molecule */
		public:
			Alkane (const int numCarbons);			// a default constructor
			virtual ~Alkane ();
			Alkane (const Molecule& molecule);		// copy constructor for casting from a molecule

			static int numAlkanes;			// total number of carbon chains in the system

			// Functions for analysis
			virtual void SetAtoms () = 0;
			void SetCarbons ();

			const Atom_it carbons_begin () const { return _carbons.begin(); }
			const Atom_it carbons_end () const { return _carbons.end(); }
			Atom_ptr_vec Carbons () const { return (_carbons); }
			AtomPtr Carbon (const int index) const { return (_carbons[index]); }

			VecR Vector_CoM_To_End ();
	};

	typedef std::vector<Alkane *> Alkane_ptr_vec;
	typedef std::vector<Alkane> Alkane_vec;

#endif
