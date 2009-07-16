#ifndef GRAPH_H_
#define GRAPH_H_

#include "utility.h"
#include "atom.h"
#include "h2o.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/property_map.hpp>


const double OHBONDLENGTH = 1.3;				// used to be 1.1
const double HBONDLENGTH  = 2.5;				// used to be 2.46
const double HBONDANGLE	= 30.0*M_PI/180.0;		// bonding angle has to be less than this value to be considered an H-bond
const double NOBONDLENGTH = 2.0;
const double NHBONDLENGTH = 1.3;		// uhmm... check this?

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

// bond types
typedef enum {unbonded, nobond, nhbond, hbond, ohbond, covalent} bondtype;


/* The idea of a connectivity matrix is the same as mapping the connections between different members of a system. In this case, we're dealing with atoms in a system, and the connections are defined by the distances between the atoms. We're going to employ a few tricks here in order to make data retrieval a bit more succinct, and also to store more information into one matrix.

Firstly, the matrix is represented by a 2D array of numbers. Each number will represent the distance between 2 atoms. The matrix is symmetric such that each row and each column represent the ordered list of atoms in the system. The upper triangle (all elements above the diagonal) hold the distances between atoms. Only using one half of a symmetric matrix ensures that we don't double our efforts and redo calculations between atoms.

Because we'll be working with water, a hydrogen bond is any bond less than 2.4 angstroms and greater than 1.3 (or so... this is to be defined below). If we want to find the number of hydrogen bonds an atom is involved in, then we look at that atom's column from the top to the diagonal, and from the diagonal to the end of the row, and count the number of bonds that fall within the distance of an H-bond.

As for the diagonal elements - instead of writing a routine that will count the number of h-bonds by looking over the rows and columns, as distances are calculated the diagonal elements for each atom will be updated to reflect the number of H-bonds formed. Thus at any time we can look at the diagonal and know immediately the number of H-bonds an atom is involved in.

This leaves the bottom-diagonal free to store more information. If two atoms are covalently bound, then the bottom diagonal element will mark this with a 1.0.
*/

using namespace std;
using namespace boost;

// Vertices are atoms
struct VertexProperties {
	Atom * atom;
	VecR position;
	std::string name;
};

// edges are bonds between atoms
struct EdgeProperties {
	EdgeProperties (const double b_length, const bondtype b_type) : distance(b_length), btype(b_type) { }
	double 		distance;
	bondtype	btype;
};

typedef boost::adjacency_list<listS, listS, undirectedS, VertexProperties, EdgeProperties> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::vertex_iterator Vertex_it;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;
typedef boost::graph_traits<Graph>::edge_iterator Edge_it;
typedef boost::graph_traits<Graph>::adjacency_iterator Adj_it;
typedef std::vector<Atom *> Atom_ptr_vec;

// generic property maps
template <class T, class Property_T> struct PropertyMap
{
	typedef typename boost::property_map<Graph, T Property_T::*>::type Type;
};

class BondGraph {

private:
	static Graph _graph;
	std::string	_sys_type;

	static PropertyMap<double,EdgeProperties>::Type 		b_length;
	static PropertyMap<bondtype,EdgeProperties>::Type 		b_type;
	static PropertyMap<Atom *,VertexProperties>::Type 		v_atom;
	static PropertyMap<VecR,VertexProperties>::Type 		v_position;
	static PropertyMap<std::string,VertexProperties>::Type	v_name;

	void _ParseAtoms (const Atom_ptr_vec& atoms);
	void _ParseBonds ();
	void _ClearBonds ();
	void _ClearAtoms ();
	void _FixSharedAtoms ();

	void _SetBond (const Vertex& vi, const Vertex& vj, const double bondlength, const bondtype btype);
	Edge _GetBond (const Vertex& vi, const Vertex& vj) const;
	Edge _GetBond (Atom const * const a1, Atom const * const a2) const;
	void _RemoveBond (const Vertex& vi, const Vertex& vj);
	void _RemoveBond (Atom const * const a1, Atom const * const a2);
	Vertex_it _FindVertex (Atom const * const ap) const;

public:

	// constructor builds the matrix based on number of atoms to analyze
	BondGraph ();
	BondGraph (const Atom_ptr_vec& atoms, std::string sys = "xyz");
	~BondGraph ();

	void SysType (std::string sys_type) { _sys_type = sys_type; }
	void UpdateGraph (const Atom_ptr_vec& atoms);

	Atom_ptr_vec BondedAtoms (
		Atom const * const ap,
		bondtype const btype = unbonded,
		std::string const name = ""
	) const;

	int NumHBonds (Atom const * const ap) const;
	int NumHBonds (Water const * const wat) const;
	coordination WaterCoordination (Water const * const wat) const;

	double Distance (Atom const * const a1, Atom const * const a2) const;

/*
	// returns the closest atoms of a given name to a given atom
	// input is the target atom's id, atomname is the name of the other atoms in the system we want returned,
	// and number is the number of nearest atoms
	// i.e. ClosestAtoms (5, O, 3) - returns the three closest O's to the atom with ID 5
	std::vector<Atom *> ClosestAtoms (const int input, const string atomname, const int number) const;
*/
};
#endif
