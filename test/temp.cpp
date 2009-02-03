//#include "../ambersystem.h"
#include "../vecr.h"
#include <boost/config.hpp>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>

using namespace boost;

class A {
public:
	double vec[3];

	A (const double x, const double y, const double z) {
		vec[0] = x;
		vec[1] = y;
		vec[2] = z;
	}

	A () {
		for (int i = 0; i < 3; i++) 
			vec[i] = 0.0;
	}
};

class V {
public:
	VecR * vec;
	std::string name;
};

class E {
public:
	double length;
	std::string type;

};
	
	typedef boost::adjacency_list<
		boost::listS, boost::listS, boost::undirectedS, 
		V, E> G;
		
	typedef G::vertex_descriptor vd;
	typedef G::edge_descriptor ed;
	typedef G::vertex_iterator v_it;
	typedef G::edge_iterator e_it;

int main () {

//	AmberSystem sys ("prmtop", "mdcrd", "mdvel");
	std::vector<VecR *> vA;
	vA.push_back (new VecR (1.0, 2.0, 3.0));
	vA.push_back (new VecR (2.0, 3.0, 4.5));
	vA.push_back (new VecR (1.5, 2.5, 6.0));

	G g (vA.size());

	// add all the vertices
	vd v1;
	v_it vi, vj, vi_end;
	
	int i = 0;
	for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
		g[*vi].vec = vA[i];
		i++;
	}

	cout << num_vertices(g) << " vertices in the graph" << endl;

	ed e;
	bool connect;
	// go through the vertex pairs to find their distances
	for (tie(vi, vi_end) = vertices(g); vi != vi_end; vi++) {
		for (vj = vi; vj != vi_end; vj++) {
			if (vj == vi) continue;

			VecR v = *g[*vi].vec - *g[*vj].vec;

			double distance = v.Magnitude();
			tie(e, connect) = add_edge(*vi, *vj, g);
			g[e].length = distance;
		}
	}

	e_it ei, ei_end;
	for (tie(ei, ei_end) = edges(g); ei != ei_end; ei++) {
		cout << g[*ei].length << endl;
	}

	
return 0;
}

/*
int main () {

	VecR size (12.0, 12.0, 20.0);
	XYZSystem sys ("pos.xyz", size, "wanniers");

	Molecule * na;
	for (int i = 0; i < 300; i++) {

		RUN(sys.Molecules()) {
			na = sys.Molecules(i);
			printf ("%s %d\n", na->Name().c_str(), na->Wanniers().size());
			//if (sys.Molecules(i)->Name() == "no3") na = sys.Molecules(i);
		}
		sys.LoadNext();
	}

return 0;
}
*/
/*
int main () {

	AmberSystem sys (PRMTOP, MDCRD, FORCE);

	vector<Atom *> int_atoms;;
	vector<Molecule *> int_mols;
	Atom * patom;
	Molecule * pmol;

	sys.LoadNext();
	RUN (sys.Molecules()) {
		pmol = sys.Molecules(i);
		double position = pmol->Atoms(0)->Y();
		if (position < 15.0) position += Atom::Size()[y];
		if (position < 70 || position > 85.0) continue;
		
		RUN2(pmol->Atoms()) {
			int_atoms.push_back (pmol->Atoms(j));
		}
		int_mols.push_back (pmol);
	}

	sys.UpdateGraph(int_atoms);
	for (int i = 300; i < 310; i++) {
		printf ("(%d) atom = ", int_atoms[i]->NumHBonds()); int_atoms[i]->Print();
		
		vector<Atom *> atoms = sys.AdjacentAtoms (int_atoms[i]);
		RUN2 (atoms) {
			printf ("%f   ", sys.Distance(int_atoms[i], atoms[j])); atoms[j]->Print();
		}
		printf ("\n");
	}

return 0;
}
*/
/*
int main () {


 	double alpha_mat [9] = {1.539, 0.0, -0.163, 0.0, 1.656, 0.0, -0.163, 0.0, 7.200};
	MatR alpha (alpha_mat);

	VecR mu (-0.058, 0.0, 0.157);

 	//double DCM_mat [9] = {0.2505, 0.0, 0.9681, 0.0, -1.0, 0.0, 0.9681, 0.0, -0.2505};
 	double DCM_mat [9] = {0.3333, 0.0, 0.9428, 0.0, -1.0, 0.0, 0.9428, 0.0, -0.3333};
	MatR DCM (DCM_mat);

	AmberSystem sys (PRMTOP, MDCRD, FORCE);
	//for (int i = 0; i < 1000; i++)
		//sys.LoadNext();

	MatR rot;
	VecR oh1;
	VecR oh2;
	VecR Z (0,0,1.0);

	RUN (sys.Molecules()) {
		if (sys.Molecules(i)->Name() != "h2o") continue;
		Water * wat = static_cast<Water *>(sys.Molecules(i));

		wat->CalcRotationData();

		std::cout << "the water is:" << std::endl;
		wat->Print();

		
		std::cout << "moving Z from 1 to lab" << std::endl;
		(wat->EulerMatrix * Z).Print();

		oh1 = *wat->OH1();
		oh2 = *wat->OH2();

		cout << "oh1" << endl;
		oh1.Print();
		(wat->EulerMatrix.Transpose() * oh1).Print();

		cout << "oh2" << endl;
		oh2.Print();
		(wat->EulerMatrix.Transpose() * oh2).Print();

		break;
	}

	return 0;
}
*/

/*
#include <iostream>
#include <boost/G/adjacency_list.hpp>

#include "../ambersystem.h"
#include "../utility.h"

#define PRMTOP	"prmtop"
#define MDCRD	"mdcrd"
#define FORCE	"mdvel"
	
using namespace boost;
	
	enum {covalent, hydrogen, none};

	struct atom {
		int ID;
	};


int main () {

	AmberSystem sys (PRMTOP, MDCRD, FORCE);

	// vd and edge properties
	typedef property<edge_weight_t, int> bond_length;		// bond length
	typedef property<edge_weight2_t, int, bond_length> bond_type;		// this is for bond type 
	typedef property<vertex_index2_t, int> atom_id;		// atom ID
	typedef property<vertex_name_t, string> atom_id;		// atom ID

	// forming the G
	typedef adjacency_list<vecS, vecS, undirectedS,
  		atom_id, bond_type> Graph;
	Graph g;

	// property objects
	property_map<Graph, vd_index2_t>::type
  		AtomID = get(vertex_index2, g);
	property_map<Graph, vd_name_t>::type
  		AtomName = get(vertex_name, g);
	property_map<Graph, edge_weight_t>::type
  		BondLength = get(edge_weight, g);
	property_map<Graph, edge_weight2_t>::type
  		BondType = get(edge_weight2, g);

	typedef G_traits<Graph>::vertex_descriptor vd;
	typedef G_traits<Graph>::edge_descriptor edge;
	typedef G_traits<Graph>::vertex_iterator v_iter;
	typedef G_traits<Graph>::edge_iterator edge_iter;
	typedef G_traits<Graph>::adjacency_iterator adj_iter;
	std::pair<vertex_iter, v_iter> vp;

	vertex atom;
	vertex_iter atom1, atom2, atoms_end;
	edge_iter ei, ei_end;
	edge bond;
	adj_iter ai, ai_end;
	
	RUN (sys) {
		atom = add_vertex(g);
		AtomID[atom] = sys[i]->ID();
		AtomName[atom] = sys[i]->Name();
	}

	for (tie(atom1, atoms_end) = vertices(g); *atom1 != 20; ++atom1) {
		for (atom2 = atom1+1; atom2 != atoms_end; ++atom2) {
			double distance = *sys[AtomID[*atom1]] - *sys[AtomID[*atom2]];

			if (distance > 5.0) continue;

			bond = add_edge(*atom1, *atom2, g).first;	
			BondLength[bond] = distance;
		}
	}

	for (tie(atom1, atoms_end) = vertices(g); *atom1 != 20; ++atom1) {
		std::cout << "atom " << AtomID[*atom1] << endl << "[";

		for (tie(ai, ai_end) = adjacent_vertices(*atom1, g); ai != ai_end; ++ai) {
			std::cout << AtomID[*ai] << "  ";
		}
		std::cout << "]" << std::endl;
	
	}



return 0;
}
*/
