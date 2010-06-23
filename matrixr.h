#ifndef MATRIXR_H_
#define MATRIXR_H_

#include <complex>
#include <iostream>
#include "vecr.h"
#include "tensor.h"

enum element {xx=0, yx=1, zx=2, xy=3, yy=4, zy=5, xz=6, yz=7, zz=8};

typedef std::vector< Double_vector >	Double_matrix;
typedef Double_matrix::iterator			Double_matrix_it;




// ***** NOTE ******
// all matrices are to be treated as column-major - see below for more info
class MatR : public tensor::tensor_t {

  public:

	MatR () : tensor::tensor_t(3,3) { }
	// constructor from a pre-built array (column-major)
	MatR (double * const elements) : tensor::tensor_t(3,3) {
	  this->Set(elements);
	}

	// A copy constructor
	MatR (const MatR& m) : tensor::tensor_t(m) { }

	~MatR () {};

	void	Set (int const row, int const col, double const val);	// Set the element
	void	Set (coord const row, coord const col, double const val);	// Set the element
	void	Set (double * const data);

	void	Set (const MatR& input);


	MatR operator* (const MatR& m) const;
	VecR operator* (const VecR& v) const;		// Vector rotation/matrix-vector inner product

};

#endif
