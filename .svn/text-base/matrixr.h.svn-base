#ifndef MATRIXR_H_
#define MATRIXR_H_

#include <complex>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>
#include "vecr.h"
#include "utility.h"

class VecR;

#undef _LINALG_

enum element {xx=0, yx, zx, xy, yy, zy, xz, yz, zz};

#ifdef _LINALG_
// pulling in C-style code from lapack to use in this c++ program

// Calculate the eigen values & vectors of a system
extern "C" { 
	int dgeev_(char const *jobvl, char const *jobvr, int const *n, double const *a,  int const *lda, double *wr, double *wi, double *vl, int const *ldvl, double *vr, int const *ldvr, double *work, int const *lwork, int *info);
}

// Calculate the inverse of a matrix
extern "C" {
	int	dgetri_ (int const *n, double *A, int const *lda, int const *ipiv, double *work, int *lwork, int *info);
}

// Calculate the L and U factorization of a matrix
extern "C" {
	int dgetrf_ (int const *m, int const *n, double *A, int const *lda, int *ipiv, int *info);
}

#endif

// ***** NOTE ******
// all matrices are to be treated as column-major - see below for more info
class MatR {

	friend class VecR;

protected:
	boost::numeric::ublas::matrix<double>	_elements;

private:
	//double _elements[9];		// elements will be entered column-major to preserve the fortran style of lapack
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
	//MatR () { _elements =boost::numeric::ublas::matrix<double> (3,3); }
	MatR () { _elements.resize(3,3); }
		//{ for (int i=0; i<9; i++) _elements[i] = 0.0; _eigenset = false; }

	// constructor from a pre-built array (column-major)
	MatR (double const * const elements) : _elements(boost::numeric::ublas::matrix<double> (3,3)) {
		for (int i = 0; i < 9; i++) {
			_elements.insert_element (i%3, i/3, elements[i]);
		}
		_eigenset = false; 
	}
	//
	// A copy constructor
	MatR (const MatR& oldMat) : _elements(boost::numeric::ublas::matrix<double> (3,3)) {
		_elements = oldMat._elements;
		_eigenset = false; 
	}

	MatR (const boost::numeric::ublas::matrix<double>& mat) : _elements(mat) {};

	~MatR () {};

// Operators
	MatR	operator+ (const MatR& input) const;		// matrix addition
	VecR 	operator* (const VecR& input) const;		// Vector rotation
	MatR 	operator* (const MatR& input) const;		// Matrix rotation
	MatR 	operator*= (const MatR& input);		// Matrix rotation
	void 	operator= (const MatR& input)		// Matrix rotation
		{ this->Set (input); }
	double	operator[] (int const index) const {	// Return the coordinate
		if (index > 8) {
			std::cout << "Trying to access illegal index in MatrixR::operator[] (int const index)\nFix This\nMatrix:" << endl;
			this->Print();
			exit(1);
		}
		return _elements(index%3, index/3); }
		
	//MatR RotateToFrame (VecR const * const frame) const;

	double	Index (int const row, int const col) const {	// Return the element
		if (row > 2 || col > 2) {
			std::cout << "Trying to access illegal index in MatrixR::Index (int const row, int const col)\nFix This\nMatrix:" << endl;
			this->Print();
			exit(1);
		}
		return _elements(row, col);
	}

	void	Set (int const row, int const col, double const val) {	// Set the element
		if (row > 2 || col > 2) {
			std::cout << "Trying to access illegal index in MatrixR::Set (int const row, int const col, double const val)\nFix This\nMatrix:" << endl;
			this->Print();
			exit(1);
		}
		{ _elements(row, col) = val; }
	}

	void	Set (double * data)
		{ for (int i=0; i<9; i++) _elements.insert_element(i%3, i/3, data[i]); }

	void	Set (const MatR& input) 
		{ _elements = input._elements; }
		
	MatR	Transpose () 	const;

#ifdef _LINALG_
	MatR 	Inverse () 	 	const;
	MatR	Diagonalize ();
#endif

// Input & matrix manipulation
	void Zero ()								// Zero all elements of a matrix
		{ for (int i=0; i<9; i++) _elements.insert_element(i%3,i/3, 0.0); }

// Output
#ifdef _LINALG_
	std::vector< complex<double> > 	EigenValues ();
	std::vector<VecR> 				EigenVectors ();
	void 						CalcEigenSystem ();
	MatR 						Quaternion ();
#endif
	double						Trace ()			const;

	void Print () const;

};

#endif
