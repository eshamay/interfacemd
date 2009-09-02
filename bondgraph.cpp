#include "bondgraph.h"

BondGraph::BondGraph () :
	_graph(0)
{ return; }

BondGraph::BondGraph (Atom_ptr_vec& atoms) :
	_graph(atoms.size())
{
	this->_ParseAtoms (atoms);
	this->_ParseDistances ();

	tie(vi, vi_end) = vertices(_graph);
	for (next = vi; vi != vi_end; vi = next) { ++next;
		_graph[*vi].position.Print();
	}


return;
}


BondGraph::~BondGraph () {

return;
}

void BondGraph::_ParseAtoms (std::vector<Atom *>& atoms) {

	// set all the vertices
	MyVertex v;
	RUN (atoms) {
		v = vertex(i, _graph);
		_graph[v].atom = atoms[i];
		_graph[v].position = atoms[i]->Position();
	}

return;
}

// Calculate the distance between each of the atoms and update the graph edge list with bonds between them
// Edges/bonds should only exist between atoms that are bound in some fashion (h-hond, covalent, etc.)
void BondGraph::_ParseDistances () {

return;
}
