#ifndef _DENSITY_H_
#define _DENSITY_H_

#include "utility.h"
#include "ambersystem.h"
#include <vector>
#include "vecr.h"
#include "atom.h"

/******** Number density calculator *********/
/* 	Given a set of atoms and start/end points on the 3 axes of a system, and the binsize, this will return a vector with the x, y and z axis number densities */
std::vector<int> NumberDensity ( AmberSystem& sys, double const start, double const end, double const binsize, const coord axis );
std::vector<int> MoleculeDensity ( AmberSystem& sys, double const start, double const end, double const binsize, const coord axis, string atomname);

#endif
