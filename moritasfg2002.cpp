#include "moritasfg2002.h"

void MoritaSFG::CalcDipole (Water * water) {

	// this plan goes along the same method of Morita/Hynes (J. Phys. Chem. B, Vol. 106, No. 3, 2002) for calculating the dipole moment. It's based on a parameterization of the water partial charges
	water->SetAtoms();

	// first we calculate the displacements of the OH bond lengths and the water angle
	double dR1 = _oh1.Magnitude() - R_eq;
	double dR2 = _oh2.Magnitude() - R_eq;
	double dTheta = acos(_oh1 < _oh2) * 180.0 / M_PI - Theta_eq;

	//printf ("%f\t%f\t%f\n", dR1, dR2, dTheta);
	//printf ("%f\t%f\t%f\n", _oh1.Magnitude(), _oh2.Magnitude(), acos(_oh1 < _oh2) * 180.0/M_PI);

	// then calculate the symmetry-adapted products
	double X1 = dR1 + dR2;
	double X2 = dR1 - dR2;
	double X3 = dTheta;

	// these are the constants from the work
	double C1 = 0.1396, C2 = -0.1196, C3 = -0.0164, C4 = -0.0288, C5 = -0.0516, C6 = 0.0532, C7 = -0.0699, C8 = 0.0169, C9 = 0.1142;

	double alpha = C7*X2 + C8*X1*X2 + C9*X2*X3;		// originally delta q1 - delta q2
	double dqO = C1*X1 + C2*X3 + C3*X1*X1 + C4*X2*X2 + C5*X3*X3 + C6*X1*X3;	// originally delta q0

	// and calculating the partial charges...
	double qO = qO_eq + dqO;		// partial charge on the oxygen
	double q1 = (alpha - qO) / 2.0;
	double q2 = (-alpha - qO) / 2.0;

	//printf ("%f\t%f\t%f\t%f\t%f\n", alpha, beta, q0, q1, q2);
	// now the dipole moment is calculated classically
	_dipole = _o->Position()*qO + _h1->Position()*q1 + _h2->Position()*q2;

return;
}


