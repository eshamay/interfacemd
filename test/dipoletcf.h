#ifndef _DIPOLETCF_H_
#define _DIPOLETCF_H_

#include <iostream>
//#include <fftw3.h>
#include "../xyzsystem.h"
#include "../utility.h"
#include "../hno3analysis.h"

using namespace std;

// define the system we're working with
//#define WATER	1
//#define	NA	1
#define slab35  1

	const string outputfilename	 = "tcf.dat";

#ifdef	WATER
	const double xSize	= 13.4724;
	const double ySize	= 15.5566;
	const double zSize	= 30.0000;
#endif

#ifdef NA
	const double xSize	= 12.4138;
	const double ySize	= 12.4138;
	const double zSize	= 12.4138;
#endif

#ifdef slab35
	const double xSize	= 12.0;
	const double ySize	= 12.0;
	const double zSize	= 20.0;
#endif

// filename of the xyzsystem
#define filename		argv[1]
#define wannierfile		argv[2]

#define OUTPUT_FREQ		50

#endif

