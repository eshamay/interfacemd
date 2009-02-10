#ifndef GRAPH_H_
#define GRAPH_H_

class Node;

/******* Edge ********/
class Edge {
	
public:
	Edge ();
	Edge (Node * i, Node * j);
	~Edge ();

	static int num_edges;

	Node * u;
	Node * v;

};

typedef std::list<Edge *>::iterator Edge_it;
typedef std::list<Edge *> Edge_ptr_list;


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

	Node * Adjacent (const Edge const * e) const;

};

typedef std::list<Node *>::iterator Node_it;
typedef std::list<Node *> Node_ptr_list;


/******* Graph ********/
class Graph {

private:
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

	int NumNodes () const { return Node::num_nodes; }
	int NumEdges () const { return Edge::num_edges; }
};

#endif
