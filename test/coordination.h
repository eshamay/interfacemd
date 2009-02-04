#ifndef	COORDINATION_H_
#define	COORDINATION_H_

#include "../utility.h"
#include "../ambersystem.h"

#define PRMTOP	"prmtop"
#define MDCRD	"mdcrd"
#define FORCE	""

#define AXIS	y

#define TIMESTEPS	10

#define POSMIN	-5.0
#define POSMAX	150.0
#define POSRES	0.25

#define INTERFACE_LOW	28.0
#define INTERFACE_HIGH	35.0
#define PBCFLIP	15.0

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

	void PrintOutput ();
	void PrintStatus ();
	void FindWaters (std::vector<Water *>& int_mols, std::vector<Atom *>& int_atoms, AmberSystem * sys);
	void FindInterfacialWaters (std::vector<Water *>& int_mols, std::vector<Atom *>& int_atoms, AmberSystem * sys);
};

#endif
