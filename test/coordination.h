#ifndef	COORDINATION_H_
#define	COORDINATION_H_

#include "../watersystem.h"

#define AVG		// define this to average the two interfaces. Make sure to define the interfaces below!
//#define RESTART
#define HIGH_COORD	OOOHHH

typedef std::vector<int> Int_vec;
typedef std::vector< vector<int> > Int_histo;
//typedef std::map<coordination, HISTOGRAM > COORD_HISTOGRAM;

class CoordinationTest : public WaterSystem {

public:

	CoordinationTest (const int argc, const char **argv, const WaterSystemParams& params);

	void OutputData ();

	void Analysis ();

protected:

	std::vector<coordination> vcoords;
	coord_map name_map;

	// For each coordination type we're going to set up a histogram along the slab to see where they mostly occur. I.e. looking at only OH-coordinated water, we'll set up the histogram to see if they are mostly found out by the interface.
	// the data structure is a map so we can easily access each histogram just by giving a coordination
	Int_histo histo;

	void InitCoordMaps ();
	void BinPosition (Water const * const wat, coordination const coord);
};


#endif
