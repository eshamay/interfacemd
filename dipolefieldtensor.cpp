#include "dipolefieldtensor.h"

DipoleFieldTensor::DipoleFieldTensor (double const * const r) {

	_r = r;		// set the pointer internally

	_distance = sqrt(_r[0]*_r[0] + _r[1]*_r[1] + _r[2]*_r[2]);	// calculate the magnitude of the distance vector
	cout << _distance << endl;

	// now we construct the tensor. Lapack may be used at some point, so we construct it in a column-major matrix format
	double C1 = 1.0/(_distance*_distance*_distance);
	double C2 = -3.0 * C1 / (_distance*_distance);

	for (int col=0; col<3; col++) {
		for (int row=0; row<3; row++) {

			_tensor[col*3 + row] = C1 + C2*_r[row]*_r[col];
		}
	}

}

DipoleFieldTensor::~DipoleFieldTensor () {

}
