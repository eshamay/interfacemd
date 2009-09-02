#ifndef _DENSITYTEST_H_
#define _DENSITYTEST_H_

#include "../watersystem.h"

using namespace std;

class AmberTest : public AmberSystem {

public:

	AmberTest (const int argc, const char **argv, const WaterSystemParams& params);

};

#endif
