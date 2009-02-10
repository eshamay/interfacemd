#include "graph.h"

// produce an empty graph
Graph::Graph () {

	_edges.clear();
	_nodes.clear();

return;
}

// remove the graph
Graph::~Graph () {
	
	Clear ();

return;
}

void Graph::Clear () {
	
	Node_it ni;
	for (ni = _nodes.begin(); ni != _nodes.end(); ni++) {
		delete (*ni);
	}
	
	Edge_it ei;
	for (ei = _edges.begin(); ei != _edges.end(); ei++) {
		delete (*ei);
	}

	_edges.clear();
	_nodes.clear();

return;
}

Edge * Graph::AddEdge (Node * u, Node * v) {
	
	Edge * e = new Edge (u, v);
	_edges.push_back (e);
	u->AddEdge (e);
	v->AddEdge (e);

return e;
}

Node * Graph::AddNode () {
	
	Node * u = new Node ();
	_nodes.push_back (u);

return u;
}

// find an edge connecting two nodes
Edge * Graph::Edge (Node * u, Node * v) {
	Edge * e = (Edge *)NULL;

	Edge_it ei;
	for (ei = _edges.begin(); ei != _edges.end(); ei++) {
		
		if ( 
			((*ei)->u == u and (*ei)->v == v) or 
			((*ei)->u == v and (*ei)->v == u) 
		   ) 
		{
				e = *ei;
		}
	}

return (e);
}

void Graph::RemoveEdge (Edge * e) {

	// first remove the edge from the graph's global edge list
	Edge_it ei;
	for (ei = _edges.begin(); ei != _edges.end(); ei++) {
		if (*ei != e) continue;

		_edges.erase(ei);
		break;
	}
	
	// and also remove it from each incident node
	Node * u = e->u;
	Node * v = e->v;
	u->RemoveEdge (e);
	v->RemoveEdge (e);

	// then remove the edge itself
	delete (e);

return;
}

// removing a node involves several steps
void Graph::RemoveNode (Node * u) {
	
	// first we remove the edges that connect the node to others
	Edge_it ei;
	for (ei = u->edges.begin(); ei != u->edges.end(); ei++) {
		RemoveEdge (*ei);
	}

	// then remove the node from the global list
	Node_it ni;
	for (ni = _nodes.begin(); ni != _nodes.end(); ni++) {
		if (*ni != u) continue;

		_nodes.erase(ni);
		break;
	}

	// then get rid of the node itself
	delete (u);

return;
}

// Find all the nodes adjacent to a given node
Node_ptr_list Graph::AdjacentNodes (const Node const * u) const {

	Node_ptr_list nodes;
	
	Edge_it ei;
	for (ei = u->edges.begin(); ei != u->edges.end(); ei++) {
		nodes.push_back (u->Adjacent(*ei));
	}

return (nodes);
}

/******* EDGES **********/

int Edge::num_edges = 0;

Edge::Edge () {

	++Edge::num_edges;
	u = (Node *)NULL;
	v = (Node *)NULL;

return;
}

Edge::Edge (Node * i, Node * j) {

	++Edge::num_edges;
	u = i;
	v = j;

return;
}

Edge::~Edge () {
	
	--Edge::num_edges;

return;
}

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

// find the node connected to this one given a particular edge
Node * Node::Adjacent (const Edge const * e) const {

	Node * u = (e->u == this) ? e->v : e->u;

return (u);
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

