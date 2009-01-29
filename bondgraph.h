#ifndef BONDGRAPH_H_
#define BONDGRAPH_H_

#include <vector>
#include <algorithm>
#include "utility.h"
#include "atom.h"
#include "h2o.h"


const double OHBONDLENGTH = 1.3;
const double HBONDLENGTH  = 2.46;
const double HBONDANGLE	= 0.866025;		// cos(theta) has to be less than this value to be considered an H-bond
const double NOBONDLENGTH = 2.0;
const double NHBONDLENGTH = 1.3;		// uhmm... check this?

// bond types
enum {covalent, hydrogen, unbonded};

// a graph is made of vertices and edges (or in our case, atoms and bonds)
class Edge;

/********** VERTEX ************/
// note that once vertices are added, they should remain static (i.e. we don't lose or add atoms to the system)
class Vertex {

public:
	Atom 			* atom;
	int				ID;
	string			name;
	std::vector<Edge *>	edges;		// all the bonds connecting to this atom

	Vertex (Atom * patom) {
		atom = patom;
		name = atom->Name();
		ID = atom->ID();
	}
};

/********** EDGE ************/
class Edge {

private:

	void SetBondType ();		// discover the bond type depending on the distance between the two atoms
	void JoinAtoms () {			// update the edge list of both atoms
		_v1->edges.push_back(this);
		_v2->edges.push_back(this);
	}

public:
	double	bondlength;		// bond length
	int		bondtype;	// bond type (covalent, hydrogen-bond, or unbonded)
	Vertex * _v1;		// the two atoms connected by this bond
	Vertex * _v2;
	
	Edge (Vertex * v1, Vertex * v2, double length) {
		_v1 = v1; _v2 = v2; bondlength = length;
		this->SetBondType ();
		this->JoinAtoms ();
	}

	~Edge ();

};


/********** BONDGRAPH ************/
class BondGraph {

private:
	std::vector<Vertex *>		_vertices;
	std::vector<Edge *>			_edges;

public:
	BondGraph ();
	BondGraph (std::vector<Atom *>& atoms);

	void UpdateGraph (std::vector<Atom *>& int_atoms);
	void ClearGraph ();

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

};

#endif
