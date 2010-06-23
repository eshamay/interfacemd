#include "matrixr.h"

/*
MatR MatR::operator+ (const MatR& input) const {
  MatR m;
  double val;
  for (int i = 0; i < 3; i++) {
	for (int j = 0; j < 3; j++) {
	  val = _matrix[i][j] + input._matrix[i][j];
	  m.Set(i,j,val);
	}}

  return (m);
}

*/

MatR MatR::operator+ (const MatR& m) const {
  MatR n(*this);
  n.plus_assign(m);
  return n;
}

MatR MatR::operator- (const MatR& m) const {
  MatR n(*this);
  n.minus_assign(m);
  return n;
}


// multiply a vector by a matrix
VecR MatR::operator* (const VecR& v) const {		// Vector rotation/matrix-vector inner product
  VecR w;
  w.assign(prod(*this, v));
  return (w);
}

// multiply a matrix by a matrix
MatR MatR::operator* (const MatR& m) const {		// Matrix rotation/multiplication
  MatR n;
  n.assign (prod(*this, m));
  return n;
}

void MatR::Set (int const row, int const col, double const val) {	// Set the element
  (*this)(row,col) = val;
}

void MatR::Set (coord const row, coord const col, double const val) {	// Set the element
  (*this)(row,col) = val;
}

void MatR::Set (const MatR& input) {
  for (unsigned i = 0; i < 3; i++) {
	for (unsigned j = 0; j < 3; j++) {
	  (*this)(i,j) = input(i,j);
	}
  }
}	// set


// set the matrix using a pre-built array of data
void MatR::Set (double * const data) {
  for (int i = 0; i < 3; i++) {
	for (int j = 0; j < 3; j++) {
	  (*this)(i,j) = data[i*3+j];
	}
  }
}


MatR MatR::Transpose() const {
  MatR n;
  n.assign(trans(*this));
  return n;
}

MatR MatR::Inverse() const {
  MatR n;
  n.assign(this->Inverse());
  return n;
}
