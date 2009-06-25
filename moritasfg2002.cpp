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


// Calculate the polarizability of a water molecule given the
MatR MoritaSFG::CalcPolarizability (Water * water) {

	// set up the internal geometry of the molecular frame, and find all the bond lengths, atoms, etc.
	MatR dcm1 (water->DCMToLabMorita(z,1));

	// this plan goes along the same method of Morita/Hynes (J. Phys. Chem. B, Vol. 106, No. 3, 2002) for calculating the dipole moment. It's based on a parameterization of the water partial atomic charges

	// first we calculate the displacements of the OH bond lengths and the water angle
	double dR1 = water->OH1()->Magnitude() - R_eq;
	double dR2 = water->OH2()->Magnitude() - R_eq;
	double dTheta = acos(*water->OH1() < *water->OH2()) * 180.0 / M_PI - Theta_eq;

	//printf ("%f\t%f\t%f\n", dR1, dR2, dTheta);
	//printf ("%f\t%f\t%f\n", _oh1.Magnitude(), _oh2.Magnitude(), acos(_oh1 < _oh2) * 180.0/M_PI);

	// first calculate the components of the polarizability of the first OH bond (in the local OH-frame)
	MatR alpha1;
	alpha1.Set(0,0, D1+D4*dR1);
	alpha1.Set(1,1, D2+D5*dR1);
	alpha1.Set(2,2, D3+D6*dR1+D7*dR1*dR1);

	MatR dcm2 (water->DCMToLabMorita(z,2));
	// now we  calculate the components of the polarizability of the second OH bond
	MatR alpha2;
	alpha2.Set(0,0, D1+D4*dR2);
	alpha2.Set(1,1, D2+D5*dR2);
	alpha2.Set(2,2, D3+D6*dR2+D7*dR2*dR2);

	// last thing to do is to rotate both the local polarizabilities into the local oh1 frame
	MatR alpha_lab = (dcm1.Transpose()*alpha1*dcm1) + (dcm2.Transpose()*alpha2*dcm2);

return (alpha_lab);
}

// update the tensor of alpha matrices for each water molecule
void MoritaSFG::UpdateAlphaTensor () {

	// first clear out the old elements and set up the big matrix
	_alpha.resize(3*_N, std::vector<double> (3*_N, 0.0));
	MatR alpha;

	for (int block_i = 0; block_i < _N; block_i++) {
	for (int block_j = 0; block_j < _N; block_j++) {
		if (block_i == block_j) {
			alpha.Set(this->CalcPolarizability(_wats[block_i]));
		}

		for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			// off-diagonal elements are zero
			if (block_i != block_j) {
				_alpha[block_i*3+i][block_j*3+j] = 0.0;
				continue;
			}

			_alpha[block_i*3+i][block_j*3+j] = alpha.Index(i,j);
		}}
	}}


/*
	for (int i = 0; i < 3*_N; i++) {
	for (int j = 0; j < 3*_N; j++) {
		printf ("% 6.2f ", _alpha[i][j]);
	}
	printf ("\n");
	}
*/

return;
}

// Calculate the dipole field tensor for the system
void MoritaSFG::UpdateDipoleFieldTensor () {

	// first clear out the dipole field tensor elements
	_T.resize(3*_N, std::vector<double> (3*_N, 0.0));

	VecR r;

	// each 3x3 block of the matrix is a sub-matrix as defined in the paper
	for (int block_i = 0; block_i < _N; block_i++) {
	for (int block_j = 0; block_j < _N; block_j++) {

		// the diagonal elements are zeros
		if (block_i == block_j)
			continue;

		// all elements are calculated from the distance vectors between the molecules
		VecR comi = _wats[block_i]->CenterOfMass();
		VecR comj = _wats[block_j]->CenterOfMass();
		r = comi.MinVector(comj,Atom::Size());

		// each sub-matrix is thus calculated from the elements of the 3x3 matrix based on distances:
		for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {

			double val = 1.0/pow(r.Magnitude(),3) - 3.0/pow(r.Magnitude(),5) * r[i]*r[j];
			_T[block_i*3+i][block_j*3+j] = val;
		}}
	}}


/*
	for (int i = 0; i < 3*_N; i++) {
	for (int j = 0; j < 3*_N; j++) {
		printf ("% 6.2f ", _T[i][j]);
	}
	printf ("\n");
	}
*/

return;
}

// calculates the effective polarizability of the system from eq. (14) in the paper
// 	alpha_eff = alpha * (1+T*alpha)^-1
// Because the _alpha and _T matrices can just be treated as 3Nx3N in size, instead of sub-matrices, let's decompose the two into large matrices for easier work.
MatR& MoritaSFG::CalcTotalPolarizability (std::vector<Water *>& wats) {

	_wats = wats;
	_N = _wats.size();
	int n = 3*_N;

	// first step is the matrix multiplication of T and alpha
	this->UpdateAlphaTensor();
	this->UpdateDipoleFieldTensor();

	// here we calculate the local field correction tensor from equation (16)
	// G = 1 + T*alpha
	// (matrix multiplication)
	double G[n*n];

	double I[n][n];		// an identity matrix
	for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {

		if (i == j)
			I[i][j] = 1.0;
		else
			I[i][j] = 0.0;

		double sum = 0.0;
		for (int k = 0; k < n; k++)
			sum += _T[i][k] * _alpha[k][j];	// T * alpha

		G[i*n+j] = I[i][j] + sum;			// 1 + T * alpha
	}}


/*
	// now construct what we need to solve equation (23)
	double h[n][3];
	for (int block_i = 0; block_i < _N; block_i++) {
		for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			h[block_i*3+i][j] = (i == j) ? 1.0 : 0.0;
		}}
	}
*/
/*
	// test the output for this funky "identity" vector
	for (int i = 0; i < n; i++) {
	for (int j = 0; j < 3; j++) {
		printf ("% 6.2f ", h[i][j]);
	}
	printf ("\n");
	}
*/

	// now calculate the inverse of 1+T*alpha to get the real G matrix.
	this->MatrixInverse(G, n);

	_A.Zero();		// the final product we're going for (equation 25)

	double val;
	for (int i = 0; i < _N; i++) {		// molecule i

		MatR alpha;						// pick out the polarizability for the given molecule i
		for (int a = 0; a < 3; a++) {
		for (int b = 0; b < 3; b++) {
			alpha.Set(a,b,_alpha[i*3+a][i*3+b]);
		}}

		for (int j = 0; j < _N; j++) {		// molecule j

			MatR g;
			for (int a = 0; a < 3; a++) {
			for (int b = 0; b < 3; b++) {
				val = G[3*i*n+a*n+3*j+b];
				g.Set(a,b,val);				// grab the specific sub matrix i,j
			}}

			_A = _A + (alpha * g);		// the sum for the particular molecule combo (i,j)
		}
	}


	//_A.Print();


/*		test that the inverse of the matrix multiplied by the matrix actually generates the identity matrix

	double J[n][n];
	// now test it out - what's H * G = G * G^-1 = ?
	for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {

		double sum = 0.0;
		for (int k = 0; k < n; k++)
			sum += G[i*n+k] * H[k][j];

		J[i][j] += sum;
	}}

	for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {
		printf ("% 6.2f ", J[i][j]);
	}
		printf ("\n");
	}
*/

return (_A);
}

// does an in-line calculation of the inverse of the square nxn matrix G
void MoritaSFG::MatrixInverse (double *G, int n) {

	int ipiv[n];
	int info;
	/*---------Call Cell LAPACK library---------*/
	dgetrf_(&n, &n, G, &n, ipiv, &info);
	if( info != 0 )
		cout << "hey! dgetrf didn't work :( - matrix factorization is dubious" << endl;

	/*---------Query workspace-------*/
	// find out how much space we need to do this thing
	double workspace;
	int tmp=-1;
	int lwork;
	dgetri_(&n, G, &n, ipiv, &workspace, &tmp, &info);

	lwork = (int)workspace;
	double work[lwork];
	//work = (double * ) malloc (sizeof(double)*lwork);
	/*---------Call Cell LAPACK library---------*/
	dgetri_(&n, G, &n, ipiv, work, &lwork, &info);
	if (info != 0)
		cout << "matrix inverse died after 2nd dgetri_ call" << endl;

	// G now holds the inverse

return;
}

VecR& MoritaSFG::CalcTotalPolarization (std::vector<Water *>& wats) {

	_p.Zero();

	RUN (wats) {
		wats[i]->CalcDipole();
		_p += wats[i]->Dipole();
	}

return _p;
}
