#include <iostream>
//#include "../utility.h"
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

/*
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
*/
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

	//sfg.AT(5,5).Print();
		//printf ("\n\nwater #%d\n", i);
		// calculate the local-mode polarizability matrix for each OH bond
		//alpha = sfg.CalcPolarizability(wat);



/*
		const VecR * oh1 = wat->OH1();
		const VecR * oh2 = wat->OH2();

		printf ("oh1 (% 8.5f) = \n", oh1->Magnitude());
		oh1->Print();
		printf ("oh2 (% 8.5f) = \n", oh2->Magnitude());
		oh2->Print();

		VecR z1 (0.0,0.0,oh1->Magnitude());
		VecR z2 (0.0,0.0,oh2->Magnitude());

		printf ("moving (0,0,1) from oh1 frame to lab\n");
		(DCM * z1).Print();
		printf ("moving oh1 from lab to local frame\n");
		(DCM.Transpose() * (*oh1)).Print();

		DCM = wat->DCMToLabMorita(z,2);
		printf ("moving (0,0,1) from oh2 frame to lab\n");
		(DCM * z2).Print();
		printf ("moving oh2 from lab to local frame\n");
		(DCM.Transpose() * (*oh2)).Print();
*/

/*
	printf ("% 12.5f% 12.5f% 12.5f\n", oh1->Magnitude(), oh2->Magnitude(), acos((*oh1) < (*oh2))*180/M_PI);
	//wat->CalcEulerAngles(y);

	//VecR oh1 = *(wat->OH1());
	VecR x (1,0,0);
	VecR y (0,1,0);
	VecR z (0,0,1);

	VecR mu (-0.058, 0.000, 0.157);
	double alpha_data[9] = {1.539, 0.000, -0.163, 0.000, 1.656, 0.000, -0.163, 0.000, 7.200};
	MatR alpha(alpha_data);
	double d[9] = {0.333, 0.0, 0.9428, 0.0, -1.0, 0.0, 0.9428, 0.0, -0.3333};
	MatR D(d);

	(dcm*alpha*dcm.Transpose()).Print();
	double sum=0.0;

	for (int i=0; i<3; i++) {
	for (int j=0; j<3; j++) {

	sum = 0.0;

	for (int l=0; l<3; l++) {
	for (int m=0; m<3; m++) {
		sum += dcm(i,l)*dcm(j,m)*alpha(l,m);
	}}

	printf ("% 12.5f", sum);
	}
	printf ("\n");
	}
	//double phi = wat->EulerAngles[0] * 180.0/M_PI;
	//double theta = wat->EulerAngles[1] * 180.0/M_PI;
	//double psi = wat->EulerAngles[2] * 180.0/M_PI;

	//printf ("% 8.3f\t% 8.3f\t% 8.3f\n", phi, theta, psi);
	//}
	//sys.LoadNext();
	//}

*/
		//printf ("% 8.3f\t% 8.3f\n", cos(2.0*atan2(1.2,1.4)), cos(2.0*atan2(1.4,1.2)));
/*
		for (int i = 0; i < 3; i++) {
			double val = wat->EulerAngles[i];
			printf ("Angle [%d] = % 5.3f\n", i, val * 180.0/M_PI);
		}

*/
	//MatR R = r_matrix (-wat->EulerAngles[0], -wat->EulerAngles[1], -wat->EulerAngles[2], 0);
	//rotate (R, wat);
	//printf ("after rotating back\n");
	//wat->Print();



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
