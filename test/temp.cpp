/*
int main () {

	AmberSystem sys ("prmtop", "mdcrd", "mdvel");
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

			g[*vi].vec->Print();
			g[*vj].vec->Print();
			if (g[*vi] == g[*vj])
				cout << "same" << endl;

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
*/
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

	AmberSystem sys ("prmtop", "mdcrd", "mdvel");

	vector<Atom *> int_atoms;;
	vector<Molecule *> int_mols;
	Atom * patom;
	Molecule * pmol;

	sys.LoadNext();
	RUN (sys.Molecules()) {
		pmol = sys.Molecules(i);
		
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

#include "../matrixr.h"
#include "../vecr.h"

int main () {

 	double da [9] = {7,9,2,8,2,4,7,1,8};
 	double db [9] = {4,3,1,9,4,7,6,4,4};
	MatR ma (da);
	MatR mb (db);


	VecR v (8,4,1);

	cout << v*v;

	cout << endl;

	(ma*mb).Print();

return 0;
}

/*
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
#include "../ambersystem.h"

#define PRMTOP	"prmtop"
#define MDCRD	"mdcrd"
#define FORCE	"mdvel"
	
int main () {

	AmberSystem sys (PRMTOP, MDCRD, FORCE);
	
	std::vector<Atom *> atoms;
	std::vector<Water *> mols;

	Molecule * mol;
	RUN (sys.Molecules()) {
		mol = sys.Molecules(i);
		if (mol->Name() != "h2o") continue;
	
		RUN2 (mol->Atoms()) {
			atoms.push_back (mol->Atoms(j));
		}
		mols.push_back (static_cast<Water *>(mol));
	}

	sys.bondgraph.UpdateGraph (atoms);

	RUN (mols) {
		coordination coord = sys.bondgraph.WaterCoordination (mols[i]);
		string name = sys.bondgraph.CoordName (coord);
		if (name == "") {
			mols[i]->Print();
		}
		printf ("%s\n", name.c_str());
	}
*/
/*
	RUN2 (atoms) {
		Atom * atom = atoms[j];
		printf ("\n");
		printf ("%s (%d)\n", atom->Name().c_str(), atom->ID());

		std::vector<Atom *> neighbors = sys.bondgraph.AdjacentAtoms (atom, ohbond);
		RUN (neighbors) {
			printf ("|\n -> %s (%d)\n", neighbors[i]->Name().c_str(), neighbors[i]->ID());
		}

		neighbors = sys.bondgraph.AdjacentAtoms (atom, hbond);
		RUN (neighbors) {
			printf ("|\n -----> %s (%d)\n", neighbors[i]->Name().c_str(), neighbors[i]->ID());
		}	
	}
return 0;
}
*/

