#include "vecr.h"


VecR::VecR() : vector_t(3) {
}


VecR::VecR (const T X, const T Y, const T Z) : vector_t(3) {
  
  (*this)(0) = X;
  (*this)(1) = Y;
  (*this)(2) = Z;
}


VecR::VecR (const VecR& oldVec) : vector_t(oldVec) {
}


VecR::VecR (const T * vec) : vector_t(3) {
  for (int i = 0; i < 3; i++)
	(*this)(i) = vec[i];
}


VecR::~VecR () {
}

VecR VecR::operator+ (const T input) const {
  VecR v;
  for (unsigned i = 0; i < 3; i++)
	v(i) = (*this)(i) + input;
  return v;
}

VecR VecR::operator+ (const VecR& input) const {
  VecR v;
  for (unsigned i = 0; i < 3; i++)
	v(i) = (*this)(i) + input(i);
  return v;
}

VecR VecR::operator- (const VecR& input) const {
  VecR v;
  for (unsigned i = 0; i < 3; i++)
	v(i) = (*this)(i) - input(i);
  return v;
}

VecR VecR::operator- (const double input) const {
  VecR v;
  for (unsigned i = 0; i < 3; i++)
	v(i) = (*this)(i) - input;
  return v;
}

void VecR::operator+= (const double input) {
  for (unsigned i = 0; i < 3; i++)
	(*this)(i) += input;
}

void VecR::operator+= (const VecR& input) {
  for (unsigned i = 0; i < 3; i++)
	(*this)(i) += input(i);
}

void VecR::operator-= (const double input) {
  return (*this -= input);
}

void VecR::operator-= (const VecR& input) {
  return (*this -= input);
}

VecR VecR::operator* (const double input) const {
  VecR v;
  for (unsigned i = 0; i < 3; i++)
	v(i) = (*this)(i) * input;
  return v;
}

// vector inner product
VecR::T VecR::operator* (const VecR& input) const {
  return inner_prod(*this, input);
}

// vector cross-product
VecR VecR::operator% (const VecR& input) const {
  VecR v;
  v(0) = (*this)(1)*input(2) - (*this)(2)*input(1);
  v(1) = (*this)(2)*input(0) - (*this)(0)*input(2);
  v(2) = (*this)(0)*input(1) - (*this)(1)*input(0);
  return (v);
}

double VecR::operator< (const VecR& input) const {
  // determine the cos(angle) between two vectors by applying the dot product, and dividing by the magnitudes
  // cos(angle) = dotproduct/magnitudes
  // Return cos(angle)
  return ( ((*this) * input) / this->Magnitude() / input.Magnitude());
}

void VecR::Zero () {
	this->clear();
}

VecR::T VecR::Magnitude () const {
  return (norm_2(*this));
}

VecR VecR::Unit () const {
  T mag = this->Magnitude();
  VecR v = *this;
  v /= mag;
  return v;
}

void VecR::Print () const {
  printf ("% 8.4e\t% 8.4e\t% 8.4e\n", (*this)(0), (*this)(1), (*this)(2));
}

