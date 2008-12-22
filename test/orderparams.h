#ifndef	ORDER_PARAMS_H_
#define	ORDER_PARAMS_H_

#include "../utility.h"
#include "../ambersystem.h"

#define PRMTOP	"prmtop"
#define MDCRD	"mdcrd"
#define FORCE	""

// set the interface correctly - different simulation = different location!!
const double INTERFACE_LOW		70.0			// the location of the interface (the top one, only, for now)
const double INTERFACE_HIGH		85.0
const double PBC_FLIP			15.0			// used for funky periodic boundaries

const coord axis = y;

const int	OUTPUT_FREQ		10					// how often the output file will be written (# of timesteps/10)
const int	TIMESTEPS		20000				// # of timesteps to process through
#define	 OUTPUTFILE		"orderparams.dat"	// name of the output file for the final spectra

// position boundaries and bin width
const double	posmin	-0.5
const double	posmax	100.0
const double	posres	1.0

const double	angmax	1.0
const double	angmin	-1.0
const double	angres	0.05

#endif
