#ifndef MATRIXR_H_
#define MATRIXR_H_

#include <complex>
#include <iostream>
#include "vecr.h"
#include "utility.h"
//#include "FTensor.h"

class VecR;

#undef _LINALG_

enum element {xx=0, yx=1, zx=2, xy=3, yy=4, zy=5, xz=6, yz=7, zz=8};

// pulling in C-style code from lapack to use in this c++ program

// Calculate the inverse of a matrix (from LAPACK)
extern "C" {
	int	dgetri_ (int const *n, double *A, int const *lda, int const *ipiv, double *work, int *lwork, int *info);
}

// Calculate the L and U factorization of a matrix
extern "C" {
	int dgetrf_ (int const *m, int const *n, double *A, int const *lda, int *ipiv, int *info);
}

// Calculate the eigen values & vectors of a system
extern "C" {
	int dgeev_(char const *jobvl, char const *jobvr, int const *n, double const *a,  int const *lda, double *wr, double *wi, double *vl, int const *ldvl, double *vr, int const *ldvr, double *work, int const *lwork, int *info);
}

//typedef FTensor::Tensor2<double,3,3>	Double_matrix;
typedef std::vector< Double_vector >	Double_matrix;

// ***** NOTE ******
// all matrices are to be treated as column-major - see below for more info
class MatR {

protected:
	Double_matrix	_matrix;

private:
	// elements will be entered column-major to preserve the fortran style of lapack when we use arrays
	/* i.e. if we were to list the indices of the elements in the matrix they would look like:
	 * 		0	3	6
	 * 		1	4	7
	 * 		2	5	8
	 */
	double _eigenvalsR[3];
	double _eigenvalsI[3];
	double _eigenvecs[9];

	bool _eigenset;

public:
	MatR () : _matrix(3, Double_vector(3)) { }

	// constructor from a pre-built array (column-major)
	MatR (double * const elements) : _matrix(3, Double_vector(3)) {
		this->Set(elements);
	}

	// constructor from a pre-built vector (column-major)
	MatR (const Double_matrix& elements) : _matrix(elements) { }

	// A copy constructor
	MatR (const MatR& oldMat) : _matrix(oldMat._matrix) { }

	~MatR () {};

// Operators
	MatR	operator+ (const MatR& input) const;		// matrix addition
	VecR 	operator* (const VecR& input) const;		// Vector rotation
	MatR 	operator* (const MatR& input) const;		// Matrix rotation
	//MatR 	operator*= (const MatR& input);				// Matrix rotation

	void 	operator= (const MatR& input) {				// assignment
		_matrix = input._matrix;
	}

	//double	operator[] (element const index) const;	// Return the coordinate
	//double	operator[] (int const index) const;	// Return the coordinate
	double	operator() (int const row, int const col) const;
	double	operator() (coord const row, coord const col) const;

	Double_matrix& Matrix() { return _matrix; }

	//MatR RotateToFrame (VecR const * const frame) const;

	double	Index (int const row, int const col) const;	// Return the element
	double	Index (coord const row, coord const col) const;	// Return the element

	void	Set (int const row, int const col, double const val);	// Set the element
	void	Set (coord const row, coord const col, double const val);	// Set the element
	void	Set (double * const data);

	void	Set (const MatR& input)
		{ _matrix = input._matrix; }

	MatR 	Inverse () 	 	const;
	MatR	Transpose () 	const;

#ifdef _LINALG_
	//MatR	Diagonalize ();
#endif

// Input & matrix manipulation
	void Zero () {								// Zero all elements of a matrix
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				_matrix[i][j] = 0.0;
			}
		}
	}

// Output
#ifdef _LINALG_
/*
	std::vector< complex<double> > 	EigenValues ();
	std::vector<VecR> 				EigenVectors ();
	void 						CalcEigenSystem ();
	MatR 						Quaternion ();
*/
#endif
	double	Trace () const;
	double Determinant () const;

	void Print () const;

};

#endif
