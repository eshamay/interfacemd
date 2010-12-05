#ifndef MOLECULE_H_
#define MOLECULE_H_

#include "vecr.h"
#include "matrixr.h"
#include "atom.h"
#include "wannier.h"
#include "utility.h"
#include <string>
#include <math.h>
#include <vector>

class Molecule {

	public:

		Molecule ();								// an empty molecule
		Molecule (const Molecule& oldMol);			// A copy constructor for performing deep copies of a molecule
		virtual ~Molecule ();

		typedef Molecule* MolPtr;
		typedef std::vector<MolPtr> Mol_ptr_vec;
		typedef Mol_ptr_vec::const_iterator Mol_it;

		typedef enum {
			NO_MOLECULE = 0,
			H, OH, H2O, H3O,
			NO3, HNO3,
			SO2,
			ALKANE, DECANE,
			CTC
		} Molecule_t;


		static int numMolecules;

		// Input functions
		void Name (std::string name) { _name = name; }	// set the molecule's name
		void MolID (int ID) { _ID = ID; }
		void MolType (Molecule::Molecule_t type) { _moltype = type; }

		// Controls
		void Shift (VecR& shift);				// Shift the origin of the entire molecule
		void clear ();							// Erases the molecule data
		VecR UpdateCenterOfMass ();				// recalculates the center of mass when coordinates are updated

		// Output Functions
		VecR CenterOfMass () const		{ return _centerofmass; }
		// A reference point within the molecule for comparing positions
		virtual VecR ReferencePoint () const { return _centerofmass; }

		/* Dealing with atoms in the molecule */
		Atom_ptr_vec Atoms () const			{ return _atoms; }

		Atom_it begin() const {
			return _atoms.begin();
		}
		Atom_it end() const {
			return _atoms.end();
		}

		const VecR_vec& Wanniers ()		const { return _wanniers; }
		VecR_it wanniers_begin ()			const { return _wanniers.begin(); }
		VecR_it wanniers_end ()				const { return _wanniers.end(); }

		double Mass ()				const 			{ return _mass; }					// Returns the molecular mass
		int size ()						const				{ return _atoms.size(); }
		std::string Name ()		const				{ return _name; }
		int MolID ()					const				{ return _ID; }
		Molecule_t MolType () const				{ return _moltype; }

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

		void Print () const;							// print out all the info of the molecule

		double MinDistance (Molecule& mol);	// calculates the minimum distance between this molecule and another (2 closest atoms)

		//VecR CalcDipole ();	// calculate the dipole
		virtual void Dipole (VecR& dip) { _dipole = dip; }
		virtual VecR Dipole () const { return _dipole; }		// return the dipole of the molecule

		virtual void Flip (const coord axis) { }
		virtual VecR MolecularAxis () { return _z; }

		// Operators
		AtomPtr operator[] (const int index) const { return _atoms[index]; }	// retrieve an atom by array index
		AtomPtr operator[] (const std::string& atomname) const;			// retrieve a particular atom using its unique name/ID
		AtomPtr operator[] (const Atom::Element_t elmt) const;
		AtomPtr GetAtom (const std::string& atomname) const;
		AtomPtr GetAtom (const Atom::Element_t elmt) const;
		//int operator+= (Atom * newAtom);					// adds an atom into the molecule

		void AddAtom (AtomPtr const newAtom);					// same as the operator

		virtual void SetAtoms () { }
		void FixAtom (AtomPtr atom);
		void FixAtoms ();
		void Rename (const std::string& name);

		MolPtr Merge (MolPtr mol);				// merges two molecules
		//int operator+= (Molecule& mol);					// Joins two molecules

		// Some stuff to work with wannier centers
		void AddWannier (const VecR& wannier) { _wanniers.push_back(wannier); } // adds a wannier center into the molecule
		void ClearWanniers () { _wanniers.clear(); }	// clear out the entire list


		// some functions to manipulate the molecule's position/orientation (symmetry operations)
		void Reflect (coord const axis, double const plane = 0.0);
		void Rotate (VecR& origin, VecR& axis, double angle);

		// get the rotation matrix to rotate a molecule to lab-frame coordinates
		virtual MatR const & DCMToLab ();
		virtual const MatR& DCM () { return _DCM; }

		template <class U>
			struct SameType : public std::binary_function<U,U,bool> {
				bool operator() (const U& left, const U& right) const {
					return left->MolType() == right->MolType();
				}
			};

		// predicate to test if the name of an atom or molecule (determined by the template parameter) is found in a vector of names
		template <class U>
			struct TypeInList : public std::binary_function<U, std::vector<Molecule_t>,bool> {
				bool operator() (const U u, const std::vector<Molecule_t>& types) const {
					return types.end() != std::find(types.begin(), types.end(), u->MolType());
				}
			};

		// returns an iterator to the first occurence of a member with the given name
		static MolPtr FindByType (const Mol_ptr_vec& mols, const Molecule_t& type) {
			Mol_it it = std::find_if (mols.begin(), mols.end(), member_functional::mem_fun_eq(&Molecule::MolType,type));
			return *it;
		}

		template <class U> 
			static void KeepByType (U& u, const Molecule_t& type) {
				u.erase(
						remove_if(u.begin(), u.end(), std::not1(member_functional::mem_fun_eq(&Molecule::MolType, type))), u.end()
						);
				return;
			}

		// keep elements of the vector with names matching one of those in the list
		template <class U> 
			static void KeepByTypes (U& u, std::vector<Molecule_t>& types) {
				u.erase(
						remove_if(u.begin(), u.end(), not1(std::bind2nd(md_name_utilities::NameInList<typename U::value_type>(), types))), u.end());
				return;
			}

		// remove all elements that have the given name
		template <class U> 
			static void RemoveByTypes (U& u, Molecule_t& type) {
				typedef const Molecule::Molecule_t& (Molecule::*fn)() const;
				u.erase(
						remove_if(u.begin(), u.end(), member_functional::mem_fun_eq(&Molecule::MolType,type)), u.end()
						);
				return;
			}

		// remove all elements that have names matching any of those in the list of names supplied
		template <class U> 
			static void RemoveByTypes (U& u, std::vector<Molecule_t>& types) {
				u.erase(
						remove_if(u.begin(), u.end(), std::bind2nd(md_name_utilities::NameInList<typename U::value_type>(), types)), u.end());
				return;
			}

		// returns an iterator to the first occurence of a member with the given name
		static MolPtr FindByID (const Mol_ptr_vec& mols, const int id) {
			Mol_it it = std::find_if (mols.begin(), mols.end(), member_functional::mem_fun_eq(&Molecule::MolID,id));
			return *it;
		}

		// predicate tells if a molecule's reference point along a given axis is above a given value
		class MoleculeAbovePosition : public std::unary_function <MolPtr,bool> {
			private:
				double position;
				coord axis;
			public:
				MoleculeAbovePosition (const double pos, const coord ax) : position(pos), axis(ax) { }
				bool operator() (const MolPtr mol) {
					return mol->ReferencePoint()[axis] > position;
				}
		};
		
		// predicate tells if a molecule's reference point along a given axis is above a given value
		class MoleculeBelowPosition : public std::unary_function <MolPtr,bool> {
			private:
				double position;
				coord axis;
			public:
				MoleculeBelowPosition (const double pos, const coord ax) : position(pos), axis(ax) { }
				bool operator() (const MolPtr mol) {
					return mol->ReferencePoint()[axis] < position;
				}
		};

	protected:
		Atom_ptr_vec	_atoms;				// the list of the atoms in the molecule
		VecR_vec		_wanniers;			// the wannier centers in the molecule
		VecR			_dipole;			// the molecular dipole
		VecR			_x, _y, _z;			// molecular frame axes


		// this is broken last I checked - not updated with coordinate updates
		VecR			_centerofmass;		// calculate by 1/M * Sum(m[i]*r[i])	where M = total mass, m[i] and r[i] are atom mass and pos

		double			_mass;				// Total molecular mass
		double			_charge;
		std::string		_name;				// some text ID or name for the molecule
		int				_ID;				// A numerical ID for the molecule
		Molecule_t		_moltype;		// an enumerated way to compare different types of molecule

		double			_eulerangles[3];	// the three euler angles theta, phi, chi
		MatR			_DCM;				// the direction cosine matrix for rotating the molecule to the lab frame
		void 			_FindEulerAngles ();// Calculates the Euler angles between the molecular axes and the fixed axes
};

typedef Molecule::MolPtr MolPtr;
typedef Molecule::Mol_ptr_vec Mol_ptr_vec;
typedef Molecule::Mol_it Mol_it;

#endif
