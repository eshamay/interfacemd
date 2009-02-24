#ifndef	COORDINATION_H_
#define	COORDINATION_H_

#include "../utility.h"
#include "../ambersystem.h"
#include "../adjacencymatrix.h"

#define AVG		// define this to average the two interfaces. Make sure to define the interfaces below!

#define PRMTOP	"prmtop"
#define MDCRD	"mdcrd"
#define FORCE	""

#define AXIS	y

#define TIMESTEPS	200000
#define RESTART		94750

#ifdef AVG
	#define POSMIN	-40.0
	#define POSMAX	40.0
#else
	#define POSMIN	-5.0
	#define POSMAX	150.0
#endif
#define POSRES	0.1

#define PBCFLIP	15.0

#define HIGH_COORD	OOOHHH

#define OUTPUT	"coord.dat"
#define OUTPUT_FREQ	25

typedef std::vector<int> HISTOGRAM;
typedef std::vector< vector<int> > COORD_HISTOGRAM;
//typedef std::map<coordination, HISTOGRAM > COORD_HISTOGRAM;

class CoordinationTest {

public:

	CoordinationTest (int argc, char **argv);

	AmberSystem sys;
	AdjacencyMatrix	matrix;

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
	void FindWaters (Water_ptr_vec& int_mols, Atom_ptr_vec& int_atoms);
	void FindInterfacialWaters (Water_ptr_vec& int_mols, Atom_ptr_vec& int_atoms);
	void BinPosition (Water const * const wat, coordination const coord);
};


#endif
