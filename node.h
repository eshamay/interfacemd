#ifndef NODE_H_
#define

#include "edge.h"

/****** Node ********/
class Node {

public:
	Node ();
	~Node ();

	Edge_ptr_list edges;		// all the edges connecting to this node

	static int num_nodes;

	void AddEdge (Edge * e);

	void RemoveEdge (Edge * e);
	void RemoveEdges ();

	int Degree () const { return edges.size(); }

};

typedef std::list<Node *>::iterator Node_it;
typedef std::list<Node *> Node_ptr_list;


#endif
