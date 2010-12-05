#include "analysis.h"

// given a system, this returns an interface location that is calculated based on the positions of either the top-most or bottom-most waters. The number of waters used in determining the interface location is given as numWaters, and the top or bottom interface is set with the bool top (true =	top surface, false = bottom)
template <typename T>
double FindAveragedWaterInterface (T& t, const int numWaters, const bool top=true);
