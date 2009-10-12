#ifndef _AMBERTEST_H_
#define _AMBERTEST_H_

#include "../watersystem.h"

using namespace std;

class AmberTest : public WaterSystem<AmberSystem> {

public:

	AmberTest (const int argc, const char **argv, const WaterSystemParams& params);

	void DoSomething ();
};

#endif
