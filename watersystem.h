#ifndef	WATERSYSTEM_H_
#define	WATERSYSTEM_H_

#include "utility.h"
#include "ambersystem.h"
#include "adjacencymatrix.h"

typedef struct {
	bool avg;

	string prmtop, mdcrd, mdvel;
	string output;

	coord axis;

	int output_freq;
	int timesteps, restart;

	double posmin, posmax, posres;
	double pbcflip;
} WaterSystemParams;

/***** DEFINITIONS *****/
/* These need to be defined in order to set several of the system parameters and data files for output */

//#define AVG		// define this to average two interfaces of a slab. Make sure to define the interfaces below!

//#define PRMTOP	// all the data files of an Amber system
//#define MDCRD
//#define FORCE		// this is modified from the MDVEL output of amber

//#define AXIS		// The primary axis that is perpendicular to the interfaces of a slab system

//#define TIMESTEPS	// number of timesteps to process
//#define RESTART	// If restarting a run - this is how many timesteps will be skipped

/* Two different position cutoffs for analysis depending on if the system is being averaged or not */
//#ifdef AVG
	//#define POSMIN
	//#define POSMAX
//#else
	//#define POSMIN
	//#define POSMAX
//#endif
//#define POSRES	// resolution of the position bins for histograms

//#define PBCFLIP	// The cutoff for the periodic boundaries. Sometimes amber screws these up and breaks up interfaces. Set this at a convenient location

//#define OUTPUT	// output file name
//#define OUTPUT_FREQ	// frequency with which something will be updated to the screen or file

class WaterSystem {

public:

	WaterSystem (const WaterSystemParams& params);
	WaterSystem (const int argc, const char **argv, const WaterSystemParams& params);
	virtual ~WaterSystem ();

	void OpenFile ();
	void OutputHeader(const WaterSystemParams& params) const;
	void Debug (string msg) const;

	void OutputStatus () const;
	void OutputData ();
	void FindWaters ();
	void SliceWaters (const double low, const double high);
	void SliceWaterCoordination (const coordination coord);
	void FindInterfacialWaters ();

	void UpdateMatrix () { matrix.UpdateMatrix (int_atoms, "sys"); }

protected:

	AmberSystem sys;
	FILE * output;

	AdjacencyMatrix	matrix;

	coord axis;

	int	output_freq;

	// position boundaries and bin width
	double	posmin, posmax, posres;
	int		posbins;
	double	pbcflip;
	unsigned int timesteps;
	unsigned int timestep, restart;

	double int_low, int_high, middle;		// the positions of analysis cutoffs


	Water_ptr_vec	int_mols;		// interfacial waters, or just all the waters in the system depending on the function call
	Atom_ptr_vec	int_atoms;		// interfacial water atoms (or as above)

};

#endif
