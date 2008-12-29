#ifndef	ORDER_PARAMS_H_
#define	ORDER_PARAMS_H_

#include "../utility.h"
#include "../ambersystem.h"

#define PRMTOP	"prmtop"
#define MDCRD	"mdcrd"
#define FORCE	""

#define INT_HIGH	74.2
#define INT_LOW		29.771

class OrderParameters {

public:

	OrderParameters (int argc, char **argv);

	AmberSystem * sys;

	FILE * output;

	vector<double> S1;
	vector<double> S2_num;
	vector<double> S2_den;
	vector<int> number_density;

	coord axis;

	int	output_freq;
	unsigned int timesteps;
	unsigned int timestep;

	// position boundaries and bin width
	double	posmin;
	double	posmax;
	double	posres;
	int		posbins;

	double  int_low, int_high, middle;

	double	angmax;
	double	angmin;
	double	angres;
	int		angbins;

	void PrintOutput ();
	void PrintStatus ();
};

#endif
