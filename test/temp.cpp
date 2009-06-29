#include <iostream>
#include "../utility.h"
//#include "../ambersystem.h"
#include "../xyzsystem.h"
#include "../moritasfg2002.h"
//#include "../pdbfile.h"

using namespace std;

/*
void movetoorigin (Molecule * wat);
MatR r_matrix (double alpha, double beta, double gamma, int fwd);
void rotate (MatR R, Molecule * mol);

*/

int main (int argc, char **argv) {
	//AmberSystem sys ("prmtop", "mdcrd", "mdvel");
	XYZSystem sys ("pos.xyz", VecR(12.0,12.0,20.0), "wanniers");
	MoritaSFG sfg;

	//for (int step = 0; step < 20; step++) {
	int steps = 25000;
	VecR M;
	std::vector<MatR> A (steps, MatR());
	MatR a;

	for (int step = 0; step < steps; step++) {
		cout << "step = " << step << endl;

		std::vector<Water *> wats;

		RUN (sys.Molecules()) {
			Molecule * mol = sys.Molecules(i);
			//if (mol->Name() != "h2o")
				//mol->Print();
			if (mol->Name() != "h2o") continue;
			Water * wat = static_cast<Water *>(mol);

			wats.push_back(wat);
		}

		if (!step)
			M = sfg.CalcTotalPolarization(wats);

		a = sfg.CalcTotalPolarizability(wats);

		A[step] = a;

		sys.LoadNext();
	}

	RUN (A) {
		double val = 0.0;
		val += A[i].Index(0,0);
		val += A[i].Index(1,1);
		val /= 2.0;

		printf ("% 15.8f% 15.8f% 15.8f\n", val*M[x], val*M[y], val*M[z]);
	}

return 0;
}
/*
int main () {

	VecR v (1.0, 2.0, 4.0);
	VecR w (2.0, 5.0, 4.0);
	VecR t = v * 0.04;
	t.Print();

	double g[9] = {
			1.0, 	2.0, 	3.0,
			4.0,	3.0,	2.0,
			8.0,	3.0,	1.0
	};
	MatR G(g);

	double h[9] = {
			4.0,	2.0,	0.0,
			1.0,	5.0,	3.0,
			4.0,	2.0,	9.0
	};
	MatR H(h);

	(G*H).Print();

return 0;
}
*/
/*

MatR r_matrix (double alpha, double beta, double gamma, int fwd) {

	// set up the rotation

	double ralpha[9] = {cos(alpha), -sin(alpha), 0.0,
					    sin(alpha), cos(alpha), 0.0,
						0.0, 			0.0, 			1.0};

	double rbeta[9] = {1.0,	0.0,		0.0,
					  0.0,	cos(beta),-sin(beta),
					  0.0,	sin(beta),cos(beta)};

	double rgamma[9] = {cos(gamma),-sin(gamma),0.0,
					  sin(gamma),cos(gamma),0.0,
					  0.0,	0.0,	1.0 };

	MatR Ralpha (ralpha);
	MatR Rbeta (rbeta);
	MatR Rgamma (rgamma);

	MatR R;
	if (fwd) {
		R = Ralpha * Rbeta * Rgamma;
	}
	else if (!fwd) {
		R = Rgamma * Rbeta * Ralpha;
	}
return R;
}

void rotate (MatR R, Molecule * mol) {
	for (int i = 0; i < 3; i++) {
		Atom * a = mol->Atoms(i);
		VecR pos = a->Position();
		pos = R * pos;
		a->Position(pos);
	}

return;
}

void movetoorigin (Molecule * wat) {

	VecR pos = wat->GetAtom("O")->Position();
	for (int i = 0; i < wat->size(); i++) {
		Atom * a = wat->Atoms(i);
		VecR atom_pos = a->Position();
		a->Position(atom_pos - pos);
	}

return;
}
*/
