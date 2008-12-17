#include "../mpi/mpisys.h"

#define PRMTOP	"prmtop"
#define MDCRD	"mdcrd"
#define FORCE	"mdvel"

using namespace std;

// set the interface correctly - different simulation = different location!!
#define INTERFACE		28.0			// the location of the interface (the top one, only, for now)

const coord axis = y;
#define	 OUTPUT_FREQ	10					// how often the output file will be written (# of timesteps/10)
#define	 TIMESTEPS		20000				// # of timesteps to process through


