#include "bondgraph.h"

using namespace boost;

/*
void Edge::SetBondType () {
	
	bondtype = unbonded;

	// Let's check OH type of bonds	
	if ( (_v1->name.find("H") != string::npos && _v2->name.find("O") != string::npos)
		|| (_v2->name.find("H") != string::npos && _v1->name.find("O") != string::npos) ) {

		if (bondlength < OHBONDLENGTH) 
			bondtype = covalent;

		// An H-bond is formed!
		if (bondlength < HBONDLENGTH && bondlength > OHBONDLENGTH) {
			bondtype = hydrogen;
			// the atoms also need to know about this so update both of them
			_v1->atom->FormHBond (_v2->atom);
			_v2->atom->FormHBond (_v1->atom);
		}
	}


	// now check for NO bonds
	if ( (_v1->name.find("N") != string::npos && _v2->name.find("O") != string::npos)
		|| (_v2->name.find("N") != string::npos && _v1->name.find("O") != string::npos) ) {

		if (bondlength < NOBONDLENGTH) 
			bondtype = covalent;
	}

	// bonding between N and H
	if ( (_v1->name.find("N") != string::npos && _v2->name.find("H") != string::npos)
		|| (_v2->name.find("N") != string::npos && _v1->name.find("H") != string::npos) ) {
	
		if (bondlength < NHBONDLENGTH) 
			bondtype = covalent;
	}



return;
}

// when removing an edge from the system, we also have to update the lists of the two vertices
Edge::~Edge () {
	
	Edge * edge = this;
	vector<Edge *>::iterator ei;

	// here we remove the this edge from the list of edges on both the vertex atoms
	ei = _v1->edges.begin();
	while (ei != _v1->edges.end()) {
		if (*ei == edge) {
			_v1->edges.erase(ei);
			break;
		}
		ei++;
	}

	ei = _v2->edges.begin();
	while (ei != _v2->edges.end()) {
		if (*ei == edge) {
			_v2->edges.erase(ei);
			break;
		}
		ei++;
	}

return;
}
*/

// default
BondGraph::BondGraph () {
	_graph.clear();

	_coord_map[UNBOUND] = "UNBOUND";
	_coord_map[O] = "O";
	_coord_map[OO] = "OO";
	_coord_map[OOO] = "OOO";
	_coord_map[H] = "H";
	_coord_map[OH] = "OH";
	_coord_map[OOH] = "OOH";
	_coord_map[OOOH] = "OOOH";
	_coord_map[HH] = "HH";
	_coord_map[OHH] = "OHH";
	_coord_map[OOHH] = "OOHH";
	_coord_map[OOOHH] = "OOOHH";
	_coord_map[HHH] = "HHH";
	_coord_map[OHHH] = "OHHH";
	_coord_map[OOHHH] = "OOHHH";
	_coord_map[OOOHHH] = "OOOHHH";

return;
}

// the constructor will set the size of the graph (i.e. number of vertices) based on the number of atoms given for analysis
BondGraph::BondGraph (std::vector<Atom *>& atoms) {
	this->UpdateGraph (atoms);
}

// runs through all the atoms and determines the new set of edges
void BondGraph::UpdateGraph (vector<Atom *>& atoms) {

	_graph.clear();

	// first let's add in all the atoms to the graph
	PATOM_IT pai, pai_end = atoms.end();
	VD vd;
	for (pai = atoms.begin(); pai != pai_end; pai++) {
		vd = add_vertex(_graph);
		_graph[vd].atom = *pai;
	}

	// first let's add in all the atoms to the graph, and also clear out their HBonding info
	RUN (atoms) {
		atoms[i]->ClearHBonds();
	}

	// now run through all atom combos and get their bondlengths
	double bondlength;
	V_IT vi, vj, vend;
	ED e;
	bool connect;
	for (tie(vi, vend) = vertices(_graph); vi != vend; vi++) {
		for (vj = vi; vj != vend; vj++) {
			if (vj == vi) continue;

			double distance = *_graph[*vi].atom - *_graph[*vj].atom;

			// check that the distance is at least an H-bond
			if (distance > HBONDLENGTH) continue;

			// add in the bond between two atoms
			tie(e, connect) = add_edge(*vi, *vj, _graph);
			_graph[e].bondlength = distance;

			// check out the bond to see what type of bond it is
			_graph[e].bond = this->BondType (*vi, *vj, e);

			// we don't, right now, have molecules where the same atom type is bonded (or H-bonded) to itself (i.e. a C-C bond, etc.)
			//if (v1->name == v2->name) continue;

			// The only bonds we're interested in are H-bonds and covalent bonds. Since covalent bonds are smaller than H-bonds, we don't need to check those.
		}
	}

return;
}

bondtype BondGraph::BondType (const VD v1, const VD v2, const ED e) const {

	bondtype bond;
	double bondlength = _graph[e].bondlength;

	// Let's check OH type of bonds	
	if ( (_graph[v1].Name("H") and _graph[v2].Name("O"))
		or (_graph[v2].Name("H") and _graph[v1].Name("O")) ) {

		if (bondlength < OHBONDLENGTH) 
			bond = ohbond;

		// An H-bond is formed!
		if (bondlength < HBONDLENGTH and bondlength > OHBONDLENGTH) {
			bond = hbond;
		}
	}

	// now check for NO bonds
	if ( (_graph[v1].Name("N") and _graph[v2].Name("O"))
		or (_graph[v2].Name("N") and _graph[v1].Name("O")) ) {

		if (bondlength < NOBONDLENGTH) 
			bond = nobond;
	}

	// bonding between N and H
	if ( (_graph[v1].Name("N") and _graph[v2].Name("O"))
		or (_graph[v2].Name("N") and _graph[v1].Name("O")) ) {
	
		if (bondlength < NHBONDLENGTH) 
			bond = nhbond;
	}

return bond;
}

V_IT BondGraph::FindVertex (const Atom * atom) const {

	// we have to find the edges associated with atom1, and see which one is connected to atom2
	V_IT vi, vend;
	for (tie(vi, vend) = vertices(_graph); vi != vend; vi++) {
		if (_graph[*vi].atom == atom)
			break;
	}

	if (vi == vend) {
		printf ("BondGraph::FindVertex () - searched for an atom in the graph that was never found. Something is wrong.\nPerhaps the list of atoms used to update the graph is not inclusive enough?\n");
		exit (1);
	}

return (vi);
}

// finds the atoms bonded (by edges) to the target atom
std::vector<Atom *> BondGraph::AdjacentAtoms (const Atom * atom) const {

	std::vector<Atom *> atoms;
	// first - find the atom in the graph
	V_IT v1 = this->FindVertex (atom);

	// then we ask the graph for all the adjacent atoms and grab them
	ADJ_IT adj, adj_end;
	for (tie(adj, adj_end) = adjacent_vertices(*v1, _graph); adj != adj_end; adj++) {
		Atom * atom = _graph[*adj].atom;
		atoms.push_back(atom);
	}

return (atoms);
}

std::vector<Atom *> BondGraph::AdjacentAtoms (const Atom * atom, const bondtype bond) const {

	std::vector<Atom *> atoms;

	// first - find the atom in the graph
	V_IT v1 = this->FindVertex (atom);
	VD v2;

	// then we ask the graph for all the bonds of the atom to analyze them to see if they match the requested bondtype
	OUT_E_IT ei, e_end;
	for (tie(ei, e_end) = out_edges(*v1, _graph); ei != e_end; ei++) {

		// check that the bond type is the one we want
		if (_graph[*ei].bond != bond) continue;
		
		v2 = target (*ei, _graph);
		atoms.push_back (_graph[v2].atom);
	}

return (atoms);
}



/*

Edge * BondGraph::FindEdge (const Atom * atom1, const Atom * atom2) const {
	
	Edge * e;
	Vertex * v1, * v2;

	// first let's find the two atoms in the graph
	v1 = this->FindVertex (atom1);
	v2 = this->FindVertex (atom2);

	// then go through all the bond combos to find the one both atoms share
	e = this->FindEdge (v1, v2);

return (e);
}

// returns the edge binding two vertices
Edge * BondGraph::FindEdge (const Vertex * v1, const Vertex * v2) const {

	Edge * e = (Edge *)NULL;
	
	// search through the edges of vertex 1 until we find the one that has v2 listed
	RUN (v1->edges) {
		RUN2 (v2->edges) {
			if (v1->edges[i] == v2->edges[j]) {
				e = v1->edges[i];
				break;
			}
		}
	}

	if (e == (Edge *)NULL) {
		printf ("BondGraph::FindEdge () hey, those two vertices aren't connected!\n Now Exiting\n");
		exit (1);
	}

return (e);
}

double BondGraph::Distance (const Atom * atom1, const Atom * atom2) const {
	Edge * e = this->FindEdge (atom1, atom2);
return (e->bondlength);
}

double BondGraph::Distance (const Vertex * atom1, const Vertex * atom2) const {
	Edge * e = this->FindEdge (atom1, atom2);
return (e->bondlength);
}

void BondGraph::ClearGraph () {
		
	// first remove all the old bonds
	vector<Edge *>::iterator ei;
	for (ei = _edges.begin(); ei != _edges.end(); ei++) {
		delete (*ei);
	}
	_edges.clear();
	
	// then remove the old atoms
	RUN(_vertices) {
		delete (_vertices[i]);
	}
	_vertices.clear();

return;
}

// find the atom at the other end of an edge from a given vertex
Atom * BondGraph::Adjacent (const Vertex * v, const Edge * e) const {
	
	Vertex * v2;

	v2 = (e->_v1 == v) ? v2 = e->_v2 : v2 = e->_v1;

return (v2->atom);
}

// finds adjacent (connected by an edge) atoms in the system with a given name
std::vector<Atom *> BondGraph::AdjacentAtoms (const Atom * atom, const string name) const {

	std::vector<Atom *> atoms (this->AdjacentAtoms (atom));
	std::vector<Atom *> matches;

	RUN (atoms) {
		if (atoms[i]->Name() != name) continue;

		matches.push_back (atoms[i]);
	}

return (matches);
}

// finds all the covalently bound atoms to a given target atom
std::vector<Atom *> BondGraph::CovalentBonds (const Atom * atom) const {

	vector<Atom *> atoms;

	// first we find the vertex
	Vertex * v = this->FindVertex (atom);

	// then grab all the atoms bonded to it that are covalent
	Edge * e;
	RUN (v->edges) {
		e = v->edges[i];
		if (e->bondtype == covalent) {
			atoms.push_back (this->Adjacent (v,e));
		}
	}

return (atoms);
}

// a more specific version of above that only looks for bound atoms of a given name
std::vector<Atom *> BondGraph::CovalentBonds (const Atom * atom, const string name) const {
	
	std::vector<Atom *> atoms (this->CovalentBonds (atom));
	std::vector<Atom *> matches;

	RUN (atoms) {
		if (atoms[i]->Name() != name) continue;

		matches.push_back (atoms[i]);
	}
			
return (matches);
}

std::vector<Atom *> BondGraph::HydrogenBonds (const Atom * atom) const {

	vector<Atom *> atoms;

	// first we find the vertex
	Vertex * v = this->FindVertex (atom);

	// then grab all the atoms bonded to it that are covalent
	Edge * e;
	RUN (v->edges) {
		e = v->edges[i];
		if (e->bondtype == hydrogen) {
			atoms.push_back (this->Adjacent (v,e));
		}
	}

return (atoms);
}

// This returns a given number of the closest atoms of a given type (name)
std::vector<Atom *> BondGraph::ClosestAtoms (const Atom * atom, const string name, const int number) const {

	std::vector<Atom *> closest;	// this is our final output
	std::vector< std::vector<double> > distances;	// used for sorting based on distance to a vertex

	Vertex * v1 (this->FindVertex (atom));
	Vertex * v2;

	// first we grab all the atoms of the given name in the system and track their distance to the target atom
	RUN (_vertices) {
		v2 = _vertices[i];
		if (v2->atom->Name().find(name) == string::npos) continue;

		// the distance between the two vertices

		std::vector<double> temp;
		temp.push_back (this->Distance (v1, v2));
		temp.push_back ((double)i);

		distances.push_back (temp);
	}
	
	// next we sort the list of the atoms we've generated based on distance
	std::sort (distances.begin(), distances.end());

	for (int i=0; i < number; i++) {
		closest.push_back (_vertices[int(distances[i][1])]->atom);
	}

return (closest);
}

coordination BondGraph::FindWaterCoordination (Water * wat) const {

	coordination coord;
	int Os = 0, Hs=0;

	// firstly we need to find all the vertices that are H-bonded to this given molecule. Thus we need to find all the atoms of the molecule.
	std::vector<Atom *> atoms = wat->Atoms();
	// now we'll check each of the atoms for H-bonding to other atoms in the system, and then count the number of bonds on the Hs and the O.
	RUN (atoms) {
		std::vector<Atom *> hbonds;
		hbonds = this->HydrogenBonds (atoms[i]);

		if (atoms[i]->Name().find("O") != string::npos) {
			Os += hbonds.size();
		}
		else if (atoms[i]->Name().find("H") != string::npos) {
			Hs += hbonds.size();
		}
	}

	coord = coordination(Os + 10*Hs);
	wat->Coordination (coord);

return coord;
}
*/
