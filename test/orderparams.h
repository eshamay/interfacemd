#ifndef	ORDER_PARAMS_H_
#define	ORDER_PARAMS_H_

#include "../utility.h"
#include "../ambersystem.h"

#define PRMTOP	"prmtop"
#define MDCRD	"mdcrd"
#define FORCE	""

class OrderParameters {

public:

	OrderParameters ();

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

	double	angmax;
	double	angmin;
	double	angres;
	int		angbins;

	void PrintOutput ();
	void PrintStatus ();
};

OrderParameters::OrderParameters () {
	
	// here is our system for analysis
	sys = new AmberSystem (PRMTOP, MDCRD, FORCE);

	output = (FILE *)NULL;
	output = fopen ("orderparams.dat", "w");
	if (output == (FILE *)NULL) {
		printf ("couldn't open the output file for reading!\n");
		exit (1);
	}

	axis = y;
	
	output_freq	=	10;					// how often the output file will be written (# of timesteps/10)
	timesteps	=	200000;				// # of timesteps to process through
	
	// position boundaries and bin width
	posmin	= -0.5;
	posmax	= 100.0;
	posres	= 0.5;
	posbins = (posmax-posmin)/posres;
	
	angmax	 = 1.0;
	angmin	 = -1.0;
	angres	 = 0.05;
	angbins = (angmax-angmin)/angres;
	
	// setup and initialize the system histogram
	// The histo looks like histo[y-position][S1][S2 numerator][S2 denominator]
	S1.clear(); S1.resize (posbins, 0.0);
	S2_num.clear(); S2_num.resize (posbins, 0.0);
	S2_den.clear(); S2_den.resize (posbins, 0.0);
	number_density.clear(); number_density.resize(posbins, 0);

return;
}
	

#endif
