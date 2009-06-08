#ifndef MORITASFG2002_H_
#define MORITASFG2002_H_

#include <complex>
#include <math.h>
//#include "ambersystem.h"
#include "h2o.h"			// defines the Water class - derived from Molecule class
#include "matrixr.h"
#include "utility.h"
#include "adjacencymatrix.h"

const double R_eq 		= 0.9575;
const double Theta_eq 	= 104.51 * M_PI/180.0;
const double Q_H_eq		= 0.3285;
const double Q_O_eq		= -0.6570;

class MoritaSFG : WaterSystem {

private:

	std:vector<VecR>	p;

	void CalcDipole (Water * water);

public:

};

#endif
