#ifndef BONDGRAPH_H_
#define BONDGRAPH_H_

#include <boost/config.hpp>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <exception>

#include "utility.h"
#include "atom.h"
#include "h2o.h"

const double OHBONDLENGTH = 1.0;
const double HBONDLENGTH  = 2.46;
const double HBONDANGLE	= 0.866025;		// cos(theta) has to be less than this value to be considered an H-bond
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

typedef std::map<coordination, string> coord_map;

// bond types
typedef enum {ohbond, nobond, nhbond, hbond, unbonded} bondtype;

/********** EDGE ************/
class EdgeDescriptor {
public:

/*
	EdgeDescriptor ();
	~EdgeDescriptor ();

	static int num_edges;
*/

	double	bondlength;	// bond length
	bondtype bond;	 	// bond type (covalent, hydrogen-bond, or unbonded)

	int id;
};

/*
int EdgeDescriptor::num_edges = 0;

EdgeDescriptor::EdgeDescriptor () {
	EdgeDescriptor::num_edges++;
	id = EdgeDescriptor::num_edges;
return;
}

EdgeDescriptor::~EdgeDescriptor () {
	EdgeDescriptor::num_edges--;
return;
}
*/


/********** VERTEX ************/
class VertexDescriptor {
public:

/*
	VertexDescriptor ();
	~VertexDescriptor ();

	static int num_vertices;
*/
	Atom * atom;

	// this can be useful to let us know if a given atom has the name we are looking for
	bool Name (std::string name) const {
		bool b;
		if (atom->Name().find(name) != std::string::npos)
			b = true;
		else
			b = false;
		return b;
	}

	int id;
};

/*
int VertexDescriptor::num_vertices = 0;

VertexDescriptor::VertexDescriptor () {
	VertexDescriptor::num_vertices++;
	id = VertexDescriptor::num_vertices;
return;
}

VertexDescriptor::~VertexDescriptor () {
	VertexDescriptor::num_vertices--;
return;
}
*/	


typedef boost::adjacency_list<
	boost::listS, boost::listS, boost::undirectedS, 
	VertexDescriptor, EdgeDescriptor> G;
		
typedef G::vertex_descriptor VD;
typedef G::edge_descriptor ED;
typedef G::vertex_iterator V_IT;
typedef G::edge_iterator E_IT;
typedef G::adjacency_iterator ADJ_IT;
typedef G::out_edge_iterator OUT_E_IT;


/********** BONDGRAPH ************/
class BondGraph {
private:
	G	_graph;		// the bondgraph object for all the data we'll need

	coord_map _coord_names;

public:

	BondGraph ();
	BondGraph (VPATOM& atoms);

	void UpdateGraph (VPATOM& int_atoms);
	void ClearGraph ();

	bondtype BondType (Atom * a1, Atom * a2, double bondlength) const;

	V_IT FindVertex (const Atom * atom) const;

	VPATOM AdjacentAtoms (const Atom * atom) const;	// finds all connected atoms (regardless of bondtype)
	VPATOM AdjacentAtoms (const Atom * atom, const bondtype bond) const;

	int NumBonds (const Atom * atom, const bondtype bond) const {
		return (AdjacentAtoms (atom, bond).size());
	}
	
	coordination WaterCoordination (const Water * wat) const;
	
	string CoordName (const coordination coord) {
		return _coord_names[coord];
	}

	coord_map& CoordNameMap () { return _coord_names; }
/*
	Edge * FindEdge (const Atom * atom1, const Atom * atom2) const;
	Edge * FindEdge (const Vertex * v1, const Vertex * v2) const;

	std::vector<Atom *> CovalentBonds (const Atom * atom) const;
	std::vector<Atom *> CovalentBonds (const Atom * atom, const string name) const;
	std::vector<Atom *> HydrogenBonds (const Atom * atom) const;

	// A routine to return the bonding coordination of a given water molecule. It will also set the coordination type in the water itself
	coordination FindWaterCoordination (Water * wat) const;

	double Distance (const Atom * atom1, const Atom * atom2) const;
	double Distance (const Vertex * atom1, const Vertex * atom2) const;

	Atom * Adjacent (const Vertex * v, const Edge * e) const;
	std::vector<Atom *> AdjacentAtoms (const Atom * atom, const string name) const;
	std::vector<Atom *> ClosestAtoms (const Atom * atom, const string name, const int number) const;
*/
};

#endif
