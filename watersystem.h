#ifndef	WATERSYSTEM_H_
#define	WATERSYSTEM_H_

#include "mdsystem.h"
#include "ambersystem.h"
#include "xyzsystem.h"
#include "utility.h"
#include "graph.h"

struct WaterSystemParams {
	WaterSystemParams (
		std::string	_output = "temp.dat",
		const int _timesteps = 200000,
		const bool _avg = false,
		const coord _axis = y,
		const VecR _ref_axis = VecR (0.0, 1.0, 0.0),
		const int _output_freq = 500, const int _restart = 0,
		const double _posmin = -20.0, const double _posmax = 150.0, 
		const double _posres = 1.0,
		const double _angmin = -1.0, const double _angmax = 1.0, 
		const double _angres = 0.01,
		const double _pbcflip = 15.0
	) :
		avg(_avg), output_filename(_output), output(fopen(_output.c_str(), "w")),
		axis(_axis), ref_axis(_ref_axis), output_freq(_output_freq), timesteps(_timesteps), restart(_restart),
		posmin(_posmin), posmax(_posmax), posres(_posres), 
		posbins ((posmax-posmin)/posres),
		pbcflip(_pbcflip),
		angmin(_angmin), angmax(_angmax), angres(_angres),
		angbins ((angmax-angmin)/angres)
	{ }

	bool avg;

	std::string output_filename;
	FILE * output;

	coord axis;
	VecR ref_axis;

	int output_freq;
	int timesteps, restart;

	double posmin, posmax, posres;
	int posbins;
	double pbcflip;
	double angmin, angmax, angres;
	int angbins;
};

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

template<class T>
class WaterSystem {

public:

	WaterSystem (const WaterSystemParams& params);

	//WaterSystem (const int argc, const char **argv, const WaterSystemParams& params);
	virtual ~WaterSystem ();

	void OpenFile ();
	void Debug (string msg) const;

	void FindWaters ();						// Find all the waters
	void FindMols (const string name);		// find all molecules with this name
	void FlipWaters (const coord axis = y);
	void SliceWaters (const double low, const double high);
	void SliceWaterCoordination (const coordination coord);
	void FindInterfacialWaters ();

	void UpdateGraph () { graph.UpdateGraph (int_atoms); }

protected:

	T * sys;
	std::string output_filename;
	FILE * output;

	BondGraph	graph;

	coord axis;
	VecR ref_axis;

	int	output_freq;

	// position boundaries and bin width
	double	posmin, posmax, posres;
	int		posbins;
	double	pbcflip;
	double 	angmin, angmax, angres;
	int		angbins;
	unsigned int timesteps;
	unsigned int restart;

	double int_low, int_high, middle;		// the positions of analysis cutoffs

	Water_ptr_vec	int_wats;		// interfacial waters, or just all the waters in the system depending on the function call
	Mol_ptr_vec int_mols;
	Atom_ptr_vec	int_atoms;		// interfacial water atoms (or as above)

};

template <class T>
WaterSystem<T>::WaterSystem (const WaterSystemParams& params)
	:
		output_filename(params.output_filename), output(params.output),
		axis(params.axis), ref_axis(params.ref_axis),
		output_freq(params.output_freq),
		posmin(params.posmin), posmax(params.posmax), posres(params.posres),
		posbins(int ((posmax - posmin)/posres) + 1), pbcflip(params.pbcflip),
		angmin(params.angmin), angmax(params.angmax), angres(params.angres),
		angbins(int ((angmax - angmin)/angres) + 1), 
		timesteps(params.timesteps), restart(params.restart)
{

/*
	if (params.avg) {	// when averaging, the interface locations are taken from the command line.
		if (argc < 3) {
			printf ("not enough parameters given (interface locations?)\n");
			exit(1);
		}
		int_low = atof(argv[1]);
		int_high = atof(argv[2]);
		middle = (int_low + int_high)/2.0;
	}
	else {
		int_low = 0.0;
		int_high = 0.0;
		middle = 0.0;
	}
*/

	OpenFile ();

	return;
}

template <class T>
WaterSystem<T>::~WaterSystem () {
	fclose (output);

	return;
}

template <class T>
void WaterSystem<T>::OpenFile () {

	if (output == (FILE *)NULL) {
		printf ("WaterSystem::WaterSystem (argc, argv) - couldn't open the data output file!\n");
		exit(1);
	}

	return;
}

template <class T>
void WaterSystem<T>::FindInterfacialWaters () {

	int_wats.clear();
	int_atoms.clear();

	Molecule * pmol;

	// go through the system
	for (int i = 0; i < sys->NumMols(); i++) {
		// grab each molecule
		pmol = sys->Molecules(i);

		// we're only looking at waters for SFG analysis right now
		if (pmol->Name() != "h2o") continue;

		// first thing we do is grab a water molecule to work with
		Water * water = static_cast<Water *>(pmol);

		double position = water->Atoms(0)->Position()[axis];
		// and find molecules that sit within the interface.
		if (position < pbcflip) position += Atom::Size()[axis];		// adjust for funky boundaries
		// these values have to be adjusted for each system
		if (position < posmin or position > posmax) continue;				// set the interface cutoffs

		int_wats.push_back (water);
		RUN2(water->Atoms()) {
			int_atoms.push_back (water->Atoms(j));
		}
	}

return;
}

template <class T>
void WaterSystem<T>::FlipWaters (const coord axis) {

	RUN (int_wats) {
		int_wats[i]->Flip(axis);
	}

return;
}

template <class T>
void WaterSystem<T>::FindMols (const string name) {

	int_mols.clear();
	int_atoms.clear();

	Molecule * pmol;

	// go through the system
	for (int i = 0; i < sys->NumMols(); i++) {
		// grab each molecule
		pmol = sys->Molecules(i);

		// we're only looking at waters for SFG analysis right now
		if (pmol->Name() != name) continue;

		// add it to the group
		int_mols.push_back (pmol);

		// and then add all of its atoms
		RUN2(pmol->Atoms()) {
			int_atoms.push_back (pmol->Atoms(j));
		}
	}

return;
}

template <class T>
void WaterSystem<T>::FindWaters () {

	int_wats.clear();
	int_atoms.clear();

	Molecule * pmol;
	Water * water;

	// go through the system
	for (int i = 0; i < sys->NumMols(); i++) {
		// grab each molecule
		pmol = sys->Molecules(i);

		// we're only looking at waters for SFG analysis right now
		if (pmol->Name() != "h2o") continue;

		// first thing we do is grab a water molecule to work with
		water = static_cast<Water *>(pmol);
		// add it to the group
		int_wats.push_back (water);

		// and then add all of its atoms
		RUN2(water->Atoms()) {
			int_atoms.push_back (water->Atoms(j));
		}
	}

return;
}


// Let's the analysis only look at a particular piece of the system instead of the entire system. That is, only use waters that are between certain positions on the long-axis
template <class T>
void WaterSystem<T>::SliceWaters (const double low, const double high) {

	std::vector<Water *> wats;
	RUN (int_wats) {
		Water * wat = int_wats[i];
		Atom * oxy = wat->GetAtom ("O");
		VecR r = oxy->Position();
		double position = r[axis];
		if (position < pbcflip) {
			position += Atom::Size()[axis];		// deal with the periodic cutoffs
		}

		if (position > low && position < high) {
			wats.push_back(wat);
		}
	}

	int_wats.clear();
	int_atoms.clear();

	RUN (wats) {
		int_wats.push_back(wats[i]);
		RUN2(wats[i]->Atoms()) {
			int_atoms.push_back(wats[i]->Atoms(j));
		}
	}

return;
}

template <class T>
void WaterSystem<T>::SliceWaterCoordination (const coordination coord) {

	Water_ptr_vec wats;

	RUN (int_wats) {
		Water * wat = int_wats[i];
		coordination c = graph.WaterCoordination(wat);
		if (c == coord) {
			wats.push_back(wat);
		}
	}

	int_wats.clear();
	int_atoms.clear();

	RUN (wats) {
		int_wats.push_back(wats[i]);
		RUN2(wats[i]->Atoms()) {
			int_atoms.push_back(wats[i]->Atoms(j));
		}
	}

return;
}

#endif
