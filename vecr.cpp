#include "vecr.h"

VecR::VecR() : _coords(3) {
}

VecR::VecR (const double X, const double Y, const double Z) : _coords(3,0.0) {
  _coords[x] = X;
  _coords[y] = Y;
  _coords[z] = Z;
}

VecR::VecR (const VecR& oldVec) : _coords(oldVec._coords) {
  //_coords = oldVec._coords;
}

VecR::VecR (const double * vec) : _coords(3,0.0) {
  for (int i = 0; i < 3; i++)
	_coords[i] = vec[i];
}

VecR::VecR (const Double_vector& vec) : _coords(vec) {
  //for (int i = 0; i < 3; i++)
  //_coords[i] = vec[i];
}

VecR::~VecR () {
}

VecR& VecR::operator= (const VecR& input) {
  _coords.resize(3,0.0);
  _coords = input._coords;
  return *this;
}

double VecR::operator[] (const coord index) const {
  return _coords[index];
}

double VecR::operator[] (const int index) const {
  if (index > 2) {
	std::cout << "Tried to access an illegal index in your vector-" << std::endl << "VecR::operator[]" << std::endl << "vector is:" << std::endl;
	this->Print();
	exit (1);
  }
  return _coords[index];
}

double VecR::operator() (const coord index) const {
  return _coords[index];
}

double VecR::operator() (const int index) const {
  return _coords[index];
}

bool VecR::operator== (const VecR& input) const {
  bool output = false;

  if (_coords[x] == input[x] && _coords[y] == input[y] && _coords[z] == input[z])
	output = true;

  return output;
}

VecR VecR::operator+ (const VecR& input) const {
  VecR v;
  for (int i = 0; i < 3; i++)
	v._coords[i] = _coords[i] + input[i];

  return (v);
}

VecR VecR::operator- (const VecR& input) const {
  VecR v;
  for (int i = 0; i < 3; i++)
	v._coords[i] = _coords[i] - input[i];
  return (v);
}

void VecR::operator+= (const VecR& input) {
  for (int i = 0; i < 3; i++)
	_coords[i] += input[i];
  return;
}

void VecR::operator+= (const double input) {
  for (int i = 0; i < 3; i++)
	_coords[i] += input;
  return;
}

// vector inner product
double VecR::operator* (const VecR& input) const {
  double val = 0.0;
  for (int i = 0; i < 3; i++)
	val += _coords[i] * input[i];
  return (val);
}

// vector scaling (multiplying by a scalar)
VecR VecR::operator* (const double input) const {
  VecR v;
  for (int i = 0; i < 3; i++)
	v._coords[i] = _coords[i] * input;
  return (v);
}

void VecR::operator*= (const double input) {
  for (int i = 0; i < 3; i++)
	_coords[i] *= input;
  return;
}

VecR VecR::operator/ (const double input) const {
  VecR v;
  for (int i = 0; i < 3; i++)
	v._coords[i] = _coords[i] / input;
  return v;
}

void VecR::operator/= (const double input) {
  for (int i = 0; i < 3; i++)
	_coords[i] /= input;
  return;
}

// vector cross-product
VecR VecR::operator% (const VecR& input) const {
  VecR v;
  v._coords[0] = _coords[1]*input._coords[2] - _coords[2]*input._coords[1];
  v._coords[1] = _coords[2]*input._coords[0] - _coords[0]*input._coords[2];
  v._coords[2] = _coords[0]*input._coords[1] - _coords[1]*input._coords[0];
  return (v);
}

void VecR::operator-= (const VecR& input) {
  for (int i = 0; i < 3; i++)
	_coords[i] -= input[i];
  return;
}

void VecR::operator-= (const double input) {
  for (int i = 0; i < 3; i++)
	_coords[i] -= input;
  return;
}

double VecR::operator< (const VecR& input) const {
  // determine the cos(angle) between two vectors by applying the dot product, and dividing by the magnitudes
  // cos(angle) = dotproduct/magnitudes
  // Return cos(angle)
  return ( ((*this) * input) / this->Magnitude() / input.Magnitude());
}

void VecR::Zero () {
  for (Double_it it = _coords.begin(); it != _coords.end(); it++) {
	*it = 0.0;
  }
}

double VecR::Magnitude () const {
  double mag = 0.0;
  for (int i=0; i<3; i++)
	mag += _coords[i]*_coords[i];
  return (sqrt(mag));
}

VecR VecR::Unit () const {
  VecR v;
  double mag = this->Magnitude();
  for (int i=0; i<3; i++)
	v._coords[i] = _coords[i] / mag;
  return (v);
}

void VecR::Print () const {
  printf ("% 8.4f\t% 8.4f\t% 8.4f\n", _coords[x], _coords[y], _coords[z]);
}

/*
// Get back a vector wrapped into a periodic cell's central-image (given by the size of the cell). this assumes that the origin is at 0,0,0
VecR& VecR::Wrap (VecR& size, VecR origin) {

while (fabs(_coords[x] - origin[x]) > size[x] / 2.0) {
if (_coords[x] < origin[x])
_coords[x] += size[x];
else
_coords[x] -= size[x];
}

while (fabs(_coords[y] - origin[y]) > size[y]/2.0) {
if (_coords[y] < origin[y])
_coords[y] += size[y];
else
_coords[y] -= size[y];
}

while (fabs(_coords[z] - origin[z]) > size[z]/2.0) {
if (_coords[z] < origin[z])
_coords[z] += size[z];
else
_coords[z] -= size[z];
}

return (*this);
}
 */
