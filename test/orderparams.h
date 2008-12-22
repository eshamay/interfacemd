#ifndef	ORDER_PARAMS_H_
#define	ORDER_PARAMS_H_

#include "../utility.h"
#include "../ambersystem.h"

#define PRMTOP	"prmtop"
#define MDCRD	"mdcrd"
#define FORCE	""

class OrderParameters {

public:

	AmberSystem * sys;

	FILE * output;

	vector< vector< vector< vector<int> > > > histo;

	vector<int> number_density;

	coord axis;

	int	output_freq;
	unsigned int timesteps;
	unsigned int timestep;

	// position boundaries and bin width
	double	posmin;
	double	posmax;
	double	posres;

	double	angmax;
	double	angmin;
	double	angres;

	void PrintOutput (int step, FILE * output, int ****histo, int *numden);
	void PrintStatus (int step);
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
	timesteps	=	20000;				// # of timesteps to process through
	output =	"orderparams.dat";		// name of the output file for the final spectra
	
	// position boundaries and bin width
	posmin	= -0.5;
	posmax	= 100.0;
	posres	= 1.0;
	const int posbins = (posmax-posmin)/posres;
	
	angmax	 = 1.0;
	angmin	 = -1.0;
	angres	 = 0.05;
	const int angbins = (angmax-angmin)/angres;
	
	// setup and initialize the system histogram
	// The histo looks like histo[y-position][S1][S2 numerator][S2 denominator]
	vector<int> a (angbins, 0);
	vector< vector<int> > b (angbins, a);
	vector< vector< vector<int> > > c (angbins, b);
	vector< vector< vector< vector<int> > > > d (posbins, c);
	histo = d;

	number_density.resize(posbins);

return;
}
	

#endif
