#ifndef GRAPH_H_
#define GRAPH_H_

#include "edge.h"
#include "node.h"

/******* Graph ********/
class Graph {

protected:
	Edge_ptr_list		_edges;
	Node_ptr_list		_nodes;

public:
	Graph ();
	~Graph ();

	Clear ();

	Edge * AddEdge (Node * u, Node * v);
	Node * AddNode ();

	void RemoveEdge (Edge * e);
	void RemoveNode (Node * e);

	Edge * Edge const (Node * u, Node * v);

	Node_ptr_list AdjacentNodes (const Node const * u) const;
	Node * Adjacent (const Node const * u, const Edge const * e) const;

	int NumNodes () const { return Node::num_nodes; }
	int NumEdges () const { return Edge::num_edges; }

};

#endif
