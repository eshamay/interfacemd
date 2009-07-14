#ifndef AMBERWATERSYSTEM_H_
#define AMBERWATERSYSTEM_H_

#include "watersystem.h"
#include "ambersystem.h"

class AmberWaterSystem : public WaterSystem<AmberSystem> {

public:
	AmberWaterSystem (const WaterSystemParams& params, std::string const prmtop, std::string const mdcrd, std::string mdvel)
	:
		WaterSystem<AmberSystem>(params)
	{
		sys = AmberSystem(prmtop, mdcrd, mdvel);
	}

protected:

	void OutputHeader () const;

};

#endif
