#include "vecr.h"

VecR::VecR () : _coords(d_vector(3,0.0)) {
}

VecR::VecR (const double X, const double Y, const double Z) : _coords(d_vector(3,0.0)) {
	_coords[x] = X;
	_coords[y] = Y;
	_coords[z] = Z;
}

VecR::VecR (const VecR& oldVec) : _coords(d_vector(3,0.0)) {
	for (int i = 0; i < 3; i++) {
		_coords[i] = oldVec._coords[i];
	}
}

VecR::VecR (const double * vec) : _coords(d_vector(3,0.0)) {
	for (int i = 0; i < 3; i++) {
		_coords[i] = vec[i];
	}
}

VecR::VecR (const d_vector& oldVec) : _coords(d_vector(3,0.0)) {
	for (int i = 0; i < 3; i++) {
		_coords[i] = oldVec[i];
	}
}

VecR::~VecR () {
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

bool VecR::operator== (const VecR& input) const {
	bool output = false;

	if (_coords[x] == input[x] && _coords[y] == input[y] && _coords[z] == input[z])
		output = true;

return output;
}

VecR VecR::operator+ (const VecR& input) const {
	VecR v (_coords);
	for (unsigned int i = 0; i < 3; i++)
		v._coords[i] += input._coords[i];
		
	return (v);
}

VecR VecR::operator- (const VecR& input) const {
	VecR v (_coords);
	for (unsigned int i = 0; i < 3; i++)
		v._coords[i] -= input._coords[i];
		
	return (v);
}

void VecR::operator+= (const VecR& input) {
	for (unsigned int i = 0; i < 3; i++)
		_coords[i] += input._coords[i];
	return;
}

double VecR::operator* (const VecR& input) const {
	double prod = 0;
	for (unsigned int i = 0; i < 3; i++) {
		prod += _coords[i] * input._coords[i];
	}
	return (prod);
}

VecR VecR::operator* (const double input) const {
	VecR v (_coords);
	for (unsigned int i = 0; i < 3; i++) {
		v._coords[i] *= input;
	}
	return (v);
}

VecR VecR::operator* (const MatR& input) const {
	VecR v;

	for (unsigned int i = 0; i < 3; i++) {
		for (unsigned int j = 0; j < 3; j++) {
			v._coords[i] += _coords[j]*input.Index(j,i);
		}
	}
	//v[x] = _coords[x]*input[xx] + _coords[y]*input[yx] + _coords[z]*input[zx];
	//v[y] = _coords[x]*input[xy] + _coords[y]*input[yy] + _coords[z]*input[zy];
	//v[z] = _coords[x]*input[xz] + _coords[y]*input[yz] + _coords[z]*input[zz];
	return (v);
}

void VecR::operator*= (const double input) {
	for (unsigned int i = 0; i < 3; i++)
		_coords[i] *= input;
	return;
}

VecR VecR::operator% (const VecR& input) const {
	VecR v;
	v._coords[x] = _coords[y]*input[z] - _coords[z]*input[y];
	v._coords[y] = _coords[z]*input[x] - _coords[x]*input[z];
	v._coords[z] = _coords[x]*input[y] - _coords[y]*input[x];

	return (v);
}

void VecR::operator-= (const VecR& input) {
	for (unsigned int i = 0; i < 3; i++)
		_coords[i] -= input[i];
	return;
}
	
double VecR::operator< (const VecR& input) const {
	// determine the cos(angle) between two vectors by applying the dot product, and dividing by the magnitudes
	// cos(angle) = dotproduct/magnitudes
	// Return cos(angle)
	return ( ((*this) * input) / this->Magnitude() / input.Magnitude());
}

void VecR::Zero () {
	for (int i=0; i<3; i++) {
		_coords[i] = 0.0;
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

// Get back a vector wrapped into a periodic cell's central-image (given by the size of the cell). this assumes that the origin is at 0,0,0
VecR VecR::Wrap (VecR size, VecR origin) {
	
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

return VecR(_coords[x], _coords[y], _coords[z]);
}

// Find the smallest vector between two locations in a periodic system defined by he size parameter
// the resulting vector will point from the current vector to the (VecR& input) location
VecR VecR::MinVector (const VecR& input, const VecR& size) const {
	
	// first we gather all our coordinates for point a (current vector) and point b (the end-point vector)
	double ax = _coords[x];
	double ay = _coords[y];
	double az = _coords[z];

	double bx = input[x];
	double by = input[y];
	double bz = input[z];
	
	double cx, cy, cz;		// this will be our output vector

	// Now we'll hold one coordinate while move the other through its periodic images until the distance between the two is <= size/2
	// this is very much like the fmod() function, but it's written out for clarity here to show that one is being held fixed
	while (fabs(ax-bx) > size[x]/2.0) {
		if (ax < bx) ax += size[x];
		else 		 ax -= size[x];
	}

	while (fabs(ay-by) > size[y]/2.0) {
		if (ay < by) ay += size[y];
		else 		 ay -= size[y];
	}

	while (fabs(az-bz) > size[z]/2.0) {
		if (az < bz) az += size[z];
		else 		 az -= size[z];
	}

	// now that both points are "closest-images" over each other, we can calculate the distance between them, and get the correct sign
	cx = bx - ax;
	cy = by - ay;
	cz = bz - az;

return (VecR(cx, cy, cz));
}

// A function for calculating the minimum-image distance between the current vector and another, given the system size
double VecR::MinDistance (const VecR& input, const VecR& size) const {
	return (MinVector(input, size).Magnitude());
}

/*
// Rotating a vector to a new coordinate axes/frame
VecR VecR::RotateToFrame (VecR const * const frame) const {
	
	VecR _x = frame[0];
	VecR _y = frame[1];
	VecR _z = frame[2];

	//here's he lab-frame coordinates hat we rotate from
	VecR X (1.0, 0.0, 0.0);
	VecR Y (0.0, 1.0, 0.0);
	VecR Z (0.0, 0.0, 1.0);

	// now we build our direction cosine matrix (eah element is he cosine of he angle between two axes)

	double rotation[9] = { _x<X, _x<Y, _x<Z, _y<X, _y<Y, _y<Z, _z<X, _z<Y, _z<Z };

	double a = rotation[0]*_coords(0) + rotation[3]*_coords(1) + rotation[6]*_coords(2);
	double b = rotation[1]*_coords(0) + rotation[4]*_coords(1) + rotation[7]*_coords(2);
	double c = rotation[2]*_coords(0) + rotation[5]*_coords(1) + rotation[8]*_coords(2);

	return (VecR(a,b,c));
}
*/
