#include "bondgraph.h"

using namespace boost;

// default
BondGraph::BondGraph () {
	_graph.clear();

	_coord_names[UNBOUND] = "UNBOUND";
	_coord_names[O] = "O";
	_coord_names[OO] = "OO";
	_coord_names[OOO] = "OOO";
	_coord_names[H] = "H";
	_coord_names[OH] = "OH";
	_coord_names[OOH] = "OOH";
	_coord_names[OOOH] = "OOOH";
	_coord_names[HH] = "HH";
	_coord_names[OHH] = "OHH";
	_coord_names[OOHH] = "OOHH";
	_coord_names[OOOHH] = "OOOHH";
	_coord_names[HHH] = "HHH";
	_coord_names[OHHH] = "OHHH";
	_coord_names[OOHHH] = "OOHHH";
	_coord_names[OOOHHH] = "OOOHHH";

return;
}

// the constructor will set the size of the graph (i.e. number of vertices) based on the number of atoms given for analysis
BondGraph::BondGraph (Atom_ptr_vec& atoms) {
	UpdateGraph (atoms);
}

// runs through all the atoms and determines the new set of edges
void BondGraph::UpdateGraph (Atom_ptr_vec& atoms) {

	Clear();

	// first let's add in all the atoms to the graph
	AtomNode * an;
	RUN (atoms) {
		an = (AtomNode *)AddNode ();
		an->atom = atoms[i];
	}

	// now run through all atom combos and get their bondlengths
	double bondlength;

	// this little for-loop bit runs through all atom-pair combinations once and find the bond-types between them
	for (unsigned int i = 0; i < _nodes.size() - 1; i++) {

		Atom * a1 = _nodes[i]->atom;

		for (unsigned int j = i + 1; j < _nodes.size(); j++) {

			Atom * a2 = _nodes[j]->atom;

			// don't connect atoms of the same name (for now.... we may later)
			if (a1->Name() == a2->Name()) continue;

			// calculate the distance between the two atoms
			double distance = *a1 - *a2;

			// check that the distance is at least an H-bond
			if (distance > HBONDLENGTH) continue;

			// add in the bond between two atoms
			AddBond (_nodes[i], _nodes[j], distance);
		}
	}

return;
}

Bond * BondGraph::AddBond (const AtomNode * u, const AtomNode * v, double length) {

	// create the new bond
	Bond * b = new Bond (u, v, length);
	// link both the atoms
	u->AddEdge (b);
	v->AddEdge (b);
	// update the bondgraph bond-list
	_edges.push_back (b);

return (b);
}

// Find a Bond that links two given atoms
Bond * BondGraph::FindBond (const Atom * x, const Atom * y) const {

	Bond * b = (Bond *)NULL;
	AtomNode * an;

	Atom * l, * m;

	for (unsigned int i = 0; i < an.size(); i++) {
		b_tmp = _edges[i];
		l = b_tmp->u->atom;
		m = b_tmp->v->atom;

		if ((l == x and m == y) or (l == y and m == x))

// unfinished
return (b);
}



V_IT BondGraph::FindAtomNode (const Atom * atom) const {

	// we have to find the edges associated with atom1, and see which one is connected to atom2
	V_IT vi, vj, vend;
	vj = (V_IT)NULL;
	for (tie(vi, vend) = vertices(_graph); vi != vend; vi++) {

		Atom * ap = _graph[*vi].atom;

		if (ap == atom)
			vj = vi;
	}

/*
	if (vj == (V_IT)NULL) {
		cout << "\nBondGraph::FindAtomNode () - searched for an atom in the graph that was never found. Something is wrong.\nPerhaps the list of atoms used to update the graph is not inclusive enough?\nThe atom is:\n" << endl;
		atom->Print();
		exit (1);
	}
*/

return (vj);
}

// finds the atoms bonded (by edges) to the target atom
Atom_ptr_vec BondGraph::AdjacentAtoms (const Atom * atom) const {

	Atom_ptr_vec atoms;
	// first - find the atom in the graph
	V_IT v1 = FindAtomNode (atom);

	// then we ask the graph for all the adjacent atoms and grab them
	ADJ_IT adj, adj_end;
	for (tie(adj, adj_end) = adjacent_vertices(*v1, _graph); adj != adj_end; adj++) {
		Atom * atom = _graph[*adj].atom;
		atoms.push_back(atom);
	}

return (atoms);
}

Atom_ptr_vec BondGraph::AdjacentAtoms (const Atom * atom, const bondtype bond) const {

	Atom_ptr_vec atoms;
	if (atom->ID() == 2490) {
		printf ("doing 2490\n");
	}

	// first - find the atom in the graph
	V_IT v1 = FindAtomNode (atom);
	if (v1 == (V_IT)NULL) {
		cout << "\nBondGraph::AdjacentAtoms - couldn't grab the atom (below) from the list given" << endl;
		atom->Print();
		printf ("Vertices in the system = %d\n", num_vertices(_graph)); fflush(stdout);
		exit(1);
	}

	// then we ask the graph for all the bonds of the atom to analyze them to see if they match the requested bondtype
	ED ed; bool connect;
	ADJ_IT adj, adj_end;
	for (tie(adj, adj_end) = adjacent_vertices(*v1, _graph); adj != adj_end; adj++) {

		// check that the bond type is the one we want
		tie(ed, connect) = edge (*v1, *adj, _graph);
		if (_graph[ed].bond != bond) continue;

//		v2 = target(*adj, _graph);

//		if (v2 == *v1) cout << "something wrong with the graph" << endl;

		atoms.push_back (_graph[*adj].atom);
	}

return (atoms);
}

coordination BondGraph::WaterCoordination (const Water * wat) const {

	if (wat->Name() != "h2o") {
		cout << "BondGraph::WaterCoordinaiton - The molecule given is named " << wat->Name() << " not \"h2o\". Is this right?" << endl;
		exit (1);
	}

	int coord = 0;

	// grab all the atoms in the water and run through them to find out the number of bonds they have
	Atom_ptr_vec atoms = wat->Atoms();
	RUN (atoms) {
		Atom * atom = atoms[i];
		// for the hydrogens add 10 for each h-bond
		if (atom->Name().find("H") != string::npos)
			coord += 10 * NumBonds(atom, hbond);
		// for oxygen add 1 for each h-bond
		else if (atom->Name().find("O") != string::npos)
			coord += NumBonds(atom, hbond);
	}

return (coordination(coord));
}

/********** BOND **********/
Bond::Bond (const AtomNode const * i, const AtomNode const * j, double length) {

	++Edge::num_edges;
	u = i;
	v = j;
	bondlength = length;

	SetBondType ();

return;
}

void Bond::SetBondType () {

	bond = unbonded;

	if (u == (AtomNode *)NULL or v == (AtomNode *)NULL) {
		cout << "BondGraph::BondType - feeding in bad atoms!\n" << endl;
		exit (1);
	}

	// Let's check for bonding between Os and Hs
	if ( (u->atom->Name().find("H") != string::npos and v->atom->Name().find("O") != string::npos)
		or (u->atom->Name().find("O") != string::npos and v->atom->Name().find("H") != string::npos) ) {

		// one type of bond is the O-H covalent
		if (bondlength < OHBONDLENGTH)
			bond = ohbond;

		// Or an H-bond is formed!
		if (bondlength < HBONDLENGTH and bondlength > OHBONDLENGTH)
			bond = hbond;
	}

return;
}

/*

Bond * BondGraph::FindBond (const Atom * atom1, const Atom * atom2) const {

	Bond * e;
	AtomNode * v1, * v2;

	// first let's find the two atoms in the graph
	v1 = this->FindAtomNode (atom1);
	v2 = this->FindAtomNode (atom2);

	// then go through all the bond combos to find the one both atoms share
	e = this->FindBond (v1, v2);

return (e);
}

// returns the edge binding two vertices
Bond * BondGraph::FindBond (const AtomNode * v1, const AtomNode * v2) const {

	Bond * e = (Bond *)NULL;

	// search through the edges of vertex 1 until we find the one that has v2 listed
	RUN (v1->edges) {
		RUN2 (v2->edges) {
			if (v1->edges[i] == v2->edges[j]) {
				e = v1->edges[i];
				break;
			}
		}
	}

	if (e == (Bond *)NULL) {
		printf ("BondGraph::FindBond () hey, those two vertices aren't connected!\n Now Exiting\n");
		exit (1);
	}

return (e);
}

double BondGraph::Distance (const Atom * atom1, const Atom * atom2) const {
	Bond * e = this->FindBond (atom1, atom2);
return (e->bondlength);
}

double BondGraph::Distance (const AtomNode * atom1, const AtomNode * atom2) const {
	Bond * e = this->FindBond (atom1, atom2);
return (e->bondlength);
}

void BondGraph::ClearGraph () {

	// first remove all the old bonds
	vector<Bond *>::iterator ei;
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
Atom * BondGraph::Adjacent (const AtomNode * v, const Bond * e) const {

	AtomNode * v2;

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
	AtomNode * v = this->FindAtomNode (atom);

	// then grab all the atoms bonded to it that are covalent
	Bond * e;
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
	AtomNode * v = this->FindAtomNode (atom);

	// then grab all the atoms bonded to it that are covalent
	Bond * e;
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

	AtomNode * v1 (this->FindAtomNode (atom));
	AtomNode * v2;

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
