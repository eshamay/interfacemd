#ifndef	COORDINATION_H_
#define	COORDINATION_H_

#include "../utility.h"
#include "../ambersystem.h"

#define PRMTOP	"prmtop"
#define MDCRD	"mdcrd"
#define FORCE	""

#define AXIS	y

#define TIMESTEPS	1000

#define POSMIN	-5.0
#define POSMAX	150.0
#define POSRES	0.1

#define INTERFACE_LOW	28.0
#define INTERFACE_HIGH	35.0
#define PBCFLIP	15.0

#define OUTPUT	"coord.dat"
#define OUTPUT_FREQ	25

typedef std::vector<int> HISTOGRAM;
typedef std::map<coordination, HISTOGRAM > COORD_HISTOGRAM;

class CoordinationTest {

public:

	CoordinationTest ();

	AmberSystem * sys;

	FILE * output;

	coord axis;

	int	output_freq;
	unsigned int timesteps;
	unsigned int timestep;

	// position boundaries and bin width
	double	posmin;
	double	posmax;
	double	posres;
	int		posbins;

	double int_low, int_high, middle;		// the positions of analysis cutoffs

	std::vector<coordination> vcoords;
	coord_map name_map;

	// For each coordination type we're going to set up a histogram along the slab to see where they mostly occur. I.e. looking at only OH-coordinated water, we'll set up the histogram to see if they are mostly found out by the interface.
	// the data structure is a map so we can easily access each histogram just by giving a coordination
	COORD_HISTOGRAM histo;

	void OutputStatus (const int step) const;
	void OutputData (const int step);
	void FindWaters (VPWATER& int_mols, VPATOM& int_atoms);
	void FindInterfacialWaters (VPWATER& int_mols, VPATOM& int_atoms);
};


#endif
