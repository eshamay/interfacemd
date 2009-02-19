#include "node.h"

/******** NODES *********/
int Node::num_nodes = 0;

// create a new unconnected node
Node::Node () {
	
	++Node::num_nodes;
	edges.clear();

return;
}

Node::~Node () {
	--Node::num_nodes;

	RemoveEdges ();

return;
}


void Node::AddEdge (Edge * e) {
	edges.push_back (e);
}

// remove a single edge from a node
void Node::RemoveEdge (Edge * e) {
	
	Edge_it ei;
	
	for (ei = edges.begin(); ei != edges.end(); ei++) {
		
		// find the matching edge in the list
		if (*ei != e) continue;

		// and then remove that edge from the list
		edges.erase(ei);
		break;
	}

return;
}

// clears out all the edges from a node
void Node::RemoveEdges () {

	edges.clear();

return;
}

