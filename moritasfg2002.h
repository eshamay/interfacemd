#ifndef MORITASFG2002_H_
#define MORITASFG2002_H_

#include <complex>
#include <math.h>
#include "utility.h"
#include "h2o.h"			// defines the Water class - derived from Molecule class

using namespace std;

const double R_eq 		= 0.9575;
const double Theta_eq 	= 104.51;
const double Q_H_eq		= 0.3285;
const double Q_O_eq		= -0.6570;

// these are the constants from the 2002 Morita work
const double C1 =  0.1396,
			 C2 = -0.1196,
			 C3 = -0.0164,
			 C4 = -0.0288,
			 C5 = -0.0516,
			 C6 =  0.0532,
			 C7 = -0.0699,
			 C8 =  0.0169,
			 C9 =  0.1142;


class MoritaSFG {

private:

	std::vector<VecR>	p;


public:

	MoritaSFG () { ; }
	VecR CalcDipole (Water * water);

};

#endif
