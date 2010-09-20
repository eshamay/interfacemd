#pragma once
#ifndef GRAPH_H_
#define GRAPH_H_

#include "mdsystem.h"

#include "atom.h"
#include "h2o.h"
#include "utility.h"

#include <map>
#include <string>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>

//#include <boost/property_map/property_map.hpp>
#include <boost/property_map.hpp>

#include <exception>

//#define ANGLE_CRITERIA

/* The idea of a connectivity matrix is the same as mapping the connections between different members of a system. In this case, we're dealing with atoms in a system, and the connections are defined by the distances between the atoms. We're going to employ a few tricks here in order to make data retrieval a bit more succinct, and also to store more information into one matrix.

	 Firstly, the matrix is represented by a 2D array of numbers. Each number will represent the distance between 2 atoms. The matrix is symmetric such that each row and each column represent the ordered list of atoms in the system. The upper triangle (all elements above the diagonal) hold the distances between atoms. Only using one half of a symmetric matrix ensures that we don't double our efforts and redo calculations between atoms.

	 Because we'll be working with water, a hydrogen bond is any bond less than 2.4 angstroms and greater than 1.3 (or so... this is to be defined below). If we want to find the number of hydrogen bonds an atom is involved in, then we look at that atom's column from the top to the diagonal, and from the diagonal to the end of the row, and count the number of bonds that fall within the distance of an H-bond.

	 As for the diagonal elements - instead of writing a routine that will count the number of h-bonds by looking over the rows and columns, as distances are calculated the diagonal elements for each atom will be updated to reflect the number of H-bonds formed. Thus at any time we can look at the diagonal and know immediately the number of H-bonds an atom is involved in.

	 This leaves the bottom-diagonal free to store more information. If two atoms are covalently bound, then the bottom diagonal element will mark this with a 1.0.
	 */

namespace bondgraph {

	using namespace boost;

	// various bondlengths to be used
	const double OHBONDLENGTH = 1.2;				// used to be 1.1
	const double HBONDLENGTH  = 2.46;				// used to be 2.46
	//const double HBONDANGLECOS	= cos(30.0*M_PI/180.0);		// bonding angle has to be bigger than this cos (i.e. smaller than ~30 degrees
	const double NOBONDLENGTH = 2.0;
	const double NHBONDLENGTH = 1.3;		// uhmm... check this?
	const double SOBONDLENGTH = 1.9;


	// bond types
	typedef enum {null = 0, unbonded, hbond, covalent} bondtype;

	// useful for tracking distances between atom pairs
	typedef std::pair<double, AtomPtr>	distance_pair;
	typedef std::vector<distance_pair> 	distance_vec;

	/* Encoding of the different coordination types
	 * The numbering is based on each O having a value of 1, and each H haveing a value of 10 (i.e. add 1 for every O, and 10 for every H...). So a water in a state of OOHH bonding would have a coordination of 22, and a coordination of 13 would be OOOH, 12 = OOH, 11 = OH, 10 = H, etc.
	 */
	typedef enum {
		UNBOUND=0, O=1, OO=2, OOO=3, OOOO=4, 			// no H
		H=10, OH=11, OOH=12, OOOH=13, OOOOH=14,			// 1 H
		HH=20, OHH=21, OOHH=22, OOOHH=23, OOOOHH=24,		// 2 Hs
		HHH=30, OHHH=31, OOHHH=32, OOOHHH=33, OOOOHHH=34,	// 3 Hs
		HHHH=40, OHHHH=41, OOHHHH=42, OOOHHHH=43, OOOOHHHH=44
	} coordination;
	// And hopefully that covers all the bonding coordination types :)

	typedef std::map<coordination, std::string> coord_map;





	class BondGraph {

		private:

			// Vertices are atoms
			struct VertexProperties {
				AtomPtr atom;
				VecR position;
				Atom::Element_t element;
			};

			// edges are bonds between atoms
			struct EdgeProperties {
				EdgeProperties (const double b_length, const bondtype b_type) : distance(b_length), btype(b_type) { }
				double 		distance;
				bondtype	btype;
			};




			typedef adjacency_list<listS, listS, undirectedS, VertexProperties, EdgeProperties> Graph;
			typedef graph_traits<Graph>::vertex_descriptor Vertex;
			typedef graph_traits<Graph>::vertex_iterator Vertex_it;
			typedef std::pair<Vertex_it, Vertex_it> Vertex_pair;
			typedef graph_traits<Graph>::edge_descriptor Edge;
			typedef graph_traits<Graph>::edge_iterator Edge_it;
			typedef graph_traits<Graph>::adjacency_iterator Adj_it;

			typedef std::vector<AtomPtr> Atom_ptr_vec;

			// generic property maps
			template <class T, class Property_T> 
				struct PropertyMap {
					typedef typename boost::property_map<Graph, T Property_T::*>::type Type;
				};


			static PropertyMap<double,EdgeProperties>::Type 			b_length;
			static PropertyMap<bondtype,EdgeProperties>::Type 		b_type;
			static PropertyMap<AtomPtr,VertexProperties>::Type 		v_atom;
			static PropertyMap<VecR,VertexProperties>::Type 			v_position;
			static PropertyMap<Atom::Element_t,VertexProperties>::Type	v_elmt;

			void _ParseAtoms (const Atom_ptr_vec& atoms);
			void _ParseBonds ();
			void _ClearBonds ();
			void _ClearAtoms ();
			void _ResolveSharedHydrogens ();

			void _SetBond (const Vertex& vi, const Vertex& vj, const double bondlength, const bondtype btype);
			Edge _GetBond (const Vertex& vi, const Vertex& vj) const;
			Edge _GetBond (const AtomPtr a1, const AtomPtr a2) const;
			void _RemoveBond (const Vertex& vi, const Vertex& vj);
			void _RemoveBond (const AtomPtr a1, const AtomPtr a2);
			Vertex_it _FindVertex (const AtomPtr atom) const;


			static Graph _graph;
			std::string _sys_type;

		public:

			// constructor builds the matrix based on number of atoms to analyze
			BondGraph ();
			BondGraph (const Atom_ptr_vec& atoms, const std::string& sys = "xyz");
			~BondGraph ();


			void SysType (std::string sys_type) { _sys_type = sys_type; }
			void UpdateGraph (const Atom_ptr_vec& atoms);

			Atom_ptr_vec BondedAtoms (
					const AtomPtr ap,
					const bondtype btype = null,
					const Atom::Element_t elmt = Atom::NO_ELEMENT
					) const;

			// Given a molecule, find the atom (of an optionally given name) that is closest to the molecule but not part of it.
			distance_pair ClosestAtom (const MolPtr&, const Atom::Element_t = Atom::NO_ELEMENT) const;
			// Given an atom, find the atom closest to it (of an optionally given name) that is not part of the same molecule
			// returns a distance pair - [distance, AtomPtr]
			distance_pair ClosestAtom (const AtomPtr, const Atom::Element_t = Atom::NO_ELEMENT, bool = false) const;

			// find the atoms that are closest to an atom
			distance_vec ClosestAtoms (const AtomPtr, const int = 1, const Atom::Element_t = Atom::NO_ELEMENT, bool = false) const;

			// find the atoms closest to a given molecule
			distance_vec ClosestAtoms (const MolPtr, const int = 1, const Atom::Element_t = Atom::NO_ELEMENT) const;

			int NumHBonds (const AtomPtr) const;
			int NumHBonds (const WaterPtr) const;
			coordination WaterCoordination (const WaterPtr) const;

			// returns the distance between two vertices in the graph
			double Distance (const Vertex&, const Vertex&) const;
			// returns the distance between two atoms
			double Distance (const AtomPtr, const AtomPtr) const;



			typedef std::exception graphex;

			struct unboundhex : public graphex {
				const char* what() const throw() { return "A hydrogen was found with no covalent bonds, and no hydrogen bonds"; }
			};

			struct multiplyboundhex : public graphex {
				const char* what() const throw() { return "A hydrogen was found with more than 1 covalent bonds"; }
			};


			class VertexIsAtom_pred : public std::binary_function<Vertex_it,AtomPtr,bool> {
				public:
					bool operator() (const Vertex_it& it, const AtomPtr& atom) const {
						return v_atom[*it] == atom;
					}
			};



			/*
			// returns the closest atoms of a given name to a given atom
			// input is the target atom's id, atomname is the name of the other atoms in the system we want returned,
			// and number is the number of nearest atoms
			// i.e. ClosestAtoms (5, O, 3) - returns the three closest O's to the atom with ID 5
			std::vector<Atom *> ClosestAtoms (const int input, const string atomname, const int number) const;
			*/
	};	// BondGraph


}	// namespace bondgraph

#endif
