#ifndef DIPOLEPARM_H_
#define DIPOLEPARM_H_

#include <stdio.h>
#include <iostream>
#include <vector>
#include "vecr.h"
#include "utility.h"


class WaterDipoleParms {

	FILE 			*_file;
	int 			_num_bins[3];		// number of bins for each of the three parameters
	double			_min[3];			// min and max values of the bins for the three parameters
	double			_max[3];
	double			_dr[3];				// bin size

	double			****_data;		// pointer arrays that will get us the data we need

	double *		_Data (double r1, double r2, double theta);

public:

	WaterDipoleParms (string parmpath);

	double Magnitude (double r1, double r2, double theta) { return _Data(r1, r2, theta)[3]; }

	VecR Dipole (double r1, double r2, double theta) {
		double x = _Data(r1, r2, theta)[0];
		double y = _Data(r1, r2, theta)[1];
		double z = _Data(r1, r2, theta)[2];
		return (VecR (x, y, z));
	}

};

#endif
