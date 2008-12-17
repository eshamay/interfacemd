#ifndef _XYZDIPOLE_ANALYSIS_H_
#define _XYZDIPOLE_ANALYSIS_H_

#include <iostream>
//#include <fftw3.h>
#include "../xyzsystem.h"
#include "../utility.h"
#include <string>
#include <sstream>

using namespace std;

// define the system we're working with
//#define WATER	1
//#define	NA	1
#define slab35  1

// define the type of run we want to do
//#define TEST_RUN	1			// a simple test that outputs the magnitude of the molecular dipoles
#define TCF_RUN	1				// Outputs the autocorrelation of the total system polarization/dipole
//#define POWER_RUN	1			// generates power spectra of individual bonds

#ifdef TEST_RUN
	const string outputfilename	 = "test.dat";
#endif

#ifdef TCF_RUN
	const string outputfilename	 = "tcf.dat";
#endif

#ifdef POWER_RUN
	const string outputfilename	 = "power.dat";
#endif

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
