#ifndef MORITASFG2002_H_
#define MORITASFG2002_H_

#include <complex>
#include <math.h>
#include "utility.h"
#include "h2o.h"			// defines the Water class - derived from Molecule class

using namespace std;

const double R_eq 		= 0.9575;		// in angstroms
const double Theta_eq 	= 104.51;		// in degrees
const double Q_H_eq		= 0.3285;		// charge units (atomic units?)
const double Q_O_eq		= -0.6570;

// these are the constants from the 2002 Morita work for calculating the dipole moment
const double C1 =  0.1396,
			 C2 = -0.1196,
			 C3 = -0.0164,
			 C4 = -0.0288,
			 C5 = -0.0516,
			 C6 =  0.0532,
			 C7 = -0.0699,
			 C8 =  0.0169,
			 C9 =  0.1142;

const double D1 = 4.6077,
			 D2 = 4.8894,
			 D3 = 5.5062,
			 D4 = 1.6890,
			 D5 = 1.6102,
			 D6 = 7.3812,
			 D7 = 3.4710;

// this is the "T" tensor found in the 2002 morita-hynes paper defined between eq's 10 and 11.
typedef std::vector< std::vector< MatR > > Ttensor;
typedef std::vector< std::vector< double > > Dbl_mtr;

class MoritaSFG {

private:

	std::vector<Water *>	_wats;	// waters of the system
	int						_N;		// number of waters
	VecR					_p;		// vector of dipoles of all the waters
	Dbl_mtr					_alpha;	// polarizabilities of the waters
	Dbl_mtr					_T;		// Dipole field tensor
	bool					_T_set;
	//double *				_G;		// Local field correction tensor
	Dbl_mtr					_alphe_eff;		// effective polarizability taking into account the local-field
	MatR					_A;

	VecR CalcDipole (Water * water);
	MatR CalcPolarizability (Water * water);

	void MatrixInverse (double *G, int n);

	void UpdateDipoleFieldTensor ();

	void UpdateAlphaTensor ();

public:

	MoritaSFG () : _T_set(false) { ; }

	MatR& CalcTotalPolarizability (std::vector<Water *>& wats);
	VecR& CalcTotalPolarization (std::vector<Water *>& wats);

};

// Calculate the inverse of a matrix (from LAPACK)
extern "C" {
	int	dgetri_ (int const *n, double *A, int const *lda, int const *ipiv, double *work, int *lwork, int *info);
}

// Calculate the L and U factorization of a matrix
extern "C" {
	int dgetrf_ (int const *m, int const *n, double *A, int const *lda, int *ipiv, int *info);
}

#endif
