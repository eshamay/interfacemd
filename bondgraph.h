#ifndef BONDGRAPH_H_
#define BONDGRAPH_H_

#include "utility.h"
#include "vecr.h"
#include "atom.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/property_map.hpp>

// Definitions of various bondlength critera and angle cutoffs
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

// bond types
typedef enum {ohbond, nobond, nhbond, hbond, unbonded, covalent} bondtype;


// vertices are atoms...
struct VertexProperties {
	Atom * atom;
	VecR position;
};

// edges are bonds between atoms
struct EdgeProperties {
	double 		distance;
	bondtype	btype;
};

using namespace boost;
using namespace std;

typedef boost::adjacency_list<listS, listS, bidirectionalS, VertexProperties, EdgeProperties> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor MyVertex;
typedef boost::graph_traits<Graph>::vertex_iterator MyVertex_it;

// here's the actual graph
class BondGraph {

private:
	Graph	_graph;

	MyVertex_it vi, vi_end, next;

	void _ParseAtoms (std::vector<Atom *>& atoms);

public:
	BondGraph ();
	BondGraph (std::vector<Atom *>& atoms);
	~BondGraph ();
};

#endif
