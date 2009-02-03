#ifndef BONDGRAPH_H_
#define BONDGRAPH_H_

#include <boost/config.hpp>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>

#include "utility.h"
#include "atom.h"
#include "h2o.h"


const double OHBONDLENGTH = 1.3;
const double HBONDLENGTH  = 2.46;
const double HBONDANGLE	= 0.866025;		// cos(theta) has to be less than this value to be considered an H-bond
const double NOBONDLENGTH = 2.0;
const double NHBONDLENGTH = 1.3;		// uhmm... check this?

// bond types
typedef enum {oh, no, nh, hbond, unbonded} bondtype;

/********** EDGE ************/
class EdgeDescriptor {
public:
	double	bondlength;	// bond length
	bondtype bond;	 	// bond type (covalent, hydrogen-bond, or unbonded)
};


/********** VERTEX ************/
// note that once vertices are added, they should remain static (i.e. we don't lose or add atoms to the system)
class VertexDescriptor {
public:
	Atom * atom;

	bool operator== (std::string name) {
		bool b;
		if (atom->Name().find(name) != std::string::npos)
			b = false;
		else
			b = true;
		return b;
	}
};

	
typedef boost::adjacency_list<
	boost::listS, boost::listS, boost::undirectedS, 
	VertexDescriptor, EdgeDescriptor> G;
		
typedef G::vertex_descriptor VD;
typedef G::edge_descriptor ED;
typedef G::vertex_iterator V_IT;
typedef G::edge_iterator E_IT;

/********** BONDGRAPH ************/
class BondGraph {
private:
	G	_graph;		// the bondgraph object for all the data we'll need

public:

	BondGraph ();
	BondGraph (std::vector<Atom *>& atoms);

	void UpdateGraph (std::vector<Atom *>& int_atoms);
	void ClearGraph () { _graph.clear(); }

/*
	Vertex * FindVertex (const Atom * atom) const;

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
	std::vector<Atom *> AdjacentAtoms (const Atom * atom) const;
	std::vector<Atom *> AdjacentAtoms (const Atom * atom, const string name) const;
	std::vector<Atom *> ClosestAtoms (const Atom * atom, const string name, const int number) const;
*/
};

#endif
