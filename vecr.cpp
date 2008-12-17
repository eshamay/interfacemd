#include "vecr.h"

VecR::VecR () : _coords(boost::numeric::ublas::vector<double> (3)){
}

VecR::VecR (const double X, const double Y, const double Z) : _coords(boost::numeric::ublas::vector<double> (3)) {
	_coords(x) = X;
	_coords(y) = Y;
	_coords(z) = Z;
}

VecR::VecR (const VecR& oldVec) : _coords(boost::numeric::ublas::vector<double> (3)) {
	_coords = boost::numeric::ublas::vector<double> (oldVec._coords);
}

VecR::VecR (const boost::numeric::ublas::vector<double>& oldVec) : _coords(boost::numeric::ublas::vector<double> (3)) { 
	_coords = oldVec;
}

VecR::~VecR () {
}

double VecR::operator[] (const coord index) const {
	return _coords(index);
}

double VecR::operator[] (const int index) const {
	if (index > 2) {
		std::cout << "Tried to access an illegal index in your vector-" << endl << "VecR::operator[]" << endl << "vector is:" << endl;
		this->Print();
		exit (1);
	}
	return _coords(index);
}

bool VecR::operator== (const VecR& input) const {
	bool output = false;

	if (_coords(x) == input[x] && _coords(y) == input[y] && _coords(z) == input[z])
		output = true;

return output;
}

VecR VecR::operator+ (const VecR& input) const {
	return VecR(_coords + input._coords);
}

VecR VecR::operator- (const VecR& input) const {
	return VecR (_coords - input._coords);
}

void VecR::operator+= (const VecR& input) {
	_coords += input._coords;
}

double VecR::operator* (const VecR& input) const {
	return inner_prod(_coords, input._coords);
}

VecR VecR::operator* (double input) const {
	return VecR(_coords * input);
}

VecR VecR::operator* (const MatR& input) const {
	return VecR (boost::numeric::ublas::prod(_coords, input._elements));
}

void VecR::operator*= (double input) {
	_coords *= input;
}

VecR VecR::operator% (const VecR& input) const {
	double i, j, k;
	i = _coords(y) * input[z] - _coords(z) * input[y];
	j = _coords(z) * input[x] - _coords(x) * input[z];
	k = _coords(x) * input[y] - _coords(y) * input[x];

	return VecR (i, j, k);
}

void VecR::operator-= (const VecR& input) {
	_coords -= input._coords;
}
	
double VecR::operator< (const VecR& input) const {
	// determine the cos(angle) between two vectors by applying the dot product, and dividing by the magnitudes
	// cos(angle) = dotproduct/magnitudes
	// Return cos(angle)

	return (inner_prod(_coords, input._coords) / norm_2(_coords) / norm_2(input._coords));

}

void VecR::Zero () {
	for (int i=0; i<3; i++) {
		_coords(i) = 0.0;
	}
}

double VecR::Magnitude () const {
	return (norm_2(_coords));
}

VecR VecR::Unit () const {
	double mag = this->Magnitude();

	return VecR (_coords(x) / mag, _coords(y) / mag, _coords(z) / mag);
}
	
void VecR::Print () const {
	printf ("% 8.4f\t% 8.4f\t% 8.4f\n", _coords(x), _coords(y), _coords(z));
}

// Get back a vector wrapped into a periodic cell's central-image (given by the size of the cell). this assumes that the origin is at 0,0,0
VecR VecR::Wrap (VecR size, VecR origin) {
	
	while (fabs(_coords(x) - origin[x]) > size[x] / 2.0) {
		if (_coords(x) < origin[x]) 	
			_coords(x) += size[x];
		else 		  			
			_coords(x) -= size[x];
	}

	while (fabs(_coords(y) - origin[y]) > size[y]/2.0) {
		if (_coords(y) < origin[y]) 	
			_coords(y) += size[y];
		else 		  			
			_coords(y) -= size[y];
	}

	while (fabs(_coords(z) - origin[z]) > size[z]/2.0) {
		if (_coords(z) < origin[z]) 	
			_coords(z) += size[z];
		else	 		  		
			_coords(z) -= size[z];
	}

return VecR(_coords(x), _coords(y), _coords(z));
/*
	if ( _coords(x) < 0.0 ) {
		_coords(x) = -1.0 * _coords(x);
		_coords(x) = fmod(_coords(x), size[x]);
		_coords(x) = size[x] - _coords(x);
	}
	else if ( _coords(x) > size[x] ) {
		_coords(x) = fmod(_coords(x), size[x]);
	}

	if ( _coords(y) < 0.0 ) {
		_coords(y) = -_coords(y);
		_coords(y) = fmod(_coords(y), size[y]);
		_coords(y) = size[y] - _coords(y);
	}
	else if ( _coords(y) > size[y] ) {
		_coords(y) = fmod(_coords(y), size[y]);
	}

	if ( _coords(z) < 0.0 ) {
		_coords(z) = -_coords(z);
		_coords(z) = fmod(_coords(z), size[z]);
		_coords(z) = size[z] - _coords(z);
	}
	else if ( _coords(z) > size[z] ) {
		_coords(z) = fmod(_coords(z), size[z]);
	}
*/
}

// Find the smallest vector between two locations in a periodic system defined by he size parameter
// the resulting vector will point from the current vector to the (VecR& input) location
VecR VecR::MinVector (const VecR& input, const VecR& size) const {
	
	// first we gather all our coordinates for point a (current vector) and point b (the end-point vector)
	double ax = _coords(x);
	double ay = _coords(y);
	double az = _coords(z);

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
	return this->MinVector(input, size).Magnitude();
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
