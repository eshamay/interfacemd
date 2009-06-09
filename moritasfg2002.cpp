#include "moritasfg2002.h"

// This returns the dipole of the water in the molecular frame.
VecR MoritaSFG::CalcDipole (Water * water) {

	water->SetAtoms();

	// this plan goes along the same method of Morita/Hynes (J. Phys. Chem. B, Vol. 106, No. 3, 2002) for calculating the dipole moment. It's based on a parameterization of the water partial atomic charges

	// first we calculate the displacements of the OH bond lengths and the water angle
	double dR1 = water->OH1()->Magnitude() - R_eq;
	double dR2 = water->OH2()->Magnitude() - R_eq;
	double dTheta = acos(*water->OH1() < *water->OH2()) * 180.0 / M_PI - Theta_eq;

	//printf ("%f\t%f\t%f\n", dR1, dR2, dTheta);
	//printf ("%f\t%f\t%f\n", _oh1.Magnitude(), _oh2.Magnitude(), acos(_oh1 < _oh2) * 180.0/M_PI);

	// then calculate the symmetry-adapted coordinates
	double X1 = dR1 + dR2;
	double X2 = dR1 - dR2;
	double X3 = dTheta;

	double dqO = C1*X1 + C2*X3 + C3*X1*X1 + C4*X2*X2 + C5*X3*X3 + C6*X1*X3;	// originally delta q0
	double f = C7*X2 + C8*X1*X2 + C9*X2*X3;		// originally delta q1 - delta q2

	// and calculating the partial atomic charges...
	double qO = Q_O_eq + dqO;		// partial charge on the oxygen
	double q1 = -Q_O_eq/2.0 - f/2.0 - dqO/2.0;	// partial charge of the 2 hydrogens
	double q2 = -Q_O_eq/2.0 + f/2.0 - dqO/2.0;

	//printf ("%f\t%f\t%f\t%f\t%f\n", f, qO, q1, q2, qO+q1+q2);
	// now the dipole moment is calculated classically, and in the OH1 molecular frame
	// 		the oxygen sits at the origin, H1 sits on the positive Z-axis,
	// 		and H2 is in the xz-plane on the positive x-side
	// 		To find H2 in the molecular frame we need the rotation matrix to move it in from the lab-fram from the lab-frame
	VecR z (0.0,0.0,1.0);
	VecR p1 = z * water->OH1()->Magnitude() * q1;		// these are in atomic units?
	VecR p2 = water->DCMToLabMorita().Transpose() * (*water->OH2()) * q2;
	VecR dipole = p1 + p2;		// this is still in the molecular OH1 frame
	// note the oxygen doesn't contribute because it's set at the origin

return dipole;
}


