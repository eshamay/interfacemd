#include <iostream>
#include "../ambersystem.h"
#include "../pdbfile.h"

using namespace std;

void movetoorigin (Molecule * wat);
MatR r_matrix (double alpha, double beta, double gamma, int fwd);
void rotate (MatR R, Molecule * mol);

int main () {

	AmberSystem sys ("prmtop", "mdcrd", "");

	// original water in molecular frame
	Water * mol = new Water();
	Atom o ("O", VecR(0.0,0.0,0.0));
	Atom h1 ("H2", VecR(0.8166,0.0,-.5771));
	Atom h2 ("H1", VecR(-0.8166,0.0,-.5771));

	mol->AddAtom(&o);
	mol->AddAtom(&h1);
	mol->AddAtom(&h2);

	vector<Molecule *> mols;

	cout << endl << "mol before rotation" << endl;
	mol->Print();

	double alpha = 90.0 * M_PI/180.0;
	double beta = 45.0 * M_PI/180.0;
	double gamma = 90.0 * M_PI/180.0;
	MatR R = r_matrix (alpha, beta, gamma, 1);
	rotate (R, mol);

	cout << endl << "mol after rotation" << endl;
	mol->Print();

	mol->SetOrderAxes();
	mol->CalcEulerAngles(z);
	for (int i = 0; i < 3; i++) {
		double val = mol->EulerAngles[i];
		printf ("Angle [%d] = % 5.3f\n", i, val * 180.0/M_PI);
	}

	/*
	Water * wat = static_cast<Water *>(sys.Molecules(180));

	Atom * w1 = wat->GetAtom("H1");
	Atom * w2 = wat->GetAtom("H2");

	movetoorigin(wat);

	printf ("before rotation\n");
	wat->Print();
	wat->SetOrderAxes();
	wat->CalcEulerAngles(z);

	for (int i = 0; i < 3; i++) {
		double val = wat->EulerAngles[i];
		printf ("Angle [%d] = % 5.3f\n", i, val * 180.0/M_PI);
	}

	MatR R = r_matrix (-wat->EulerAngles[0], -wat->EulerAngles[1], -wat->EulerAngles[2], 0);
	rotate (R, wat);
	printf ("after rotating back\n");
	wat->Print();
*/



return 0;
}


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
