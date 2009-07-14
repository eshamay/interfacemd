#ifndef XYZWATERSYSTEM_H_
#define XYZWATERSYSTEM_H_

#include "watersystem.h"
#include "xyzsystem.h"

class XYZWaterSystem : public WaterSystem<XYZSystem> {

public:
	XYZWaterSystem (const WaterSystemParams& params, std::string const xyzfilepath, const VecR& size, std::string const wanniers = "")
	:
		WaterSystem<XYZSystem>(params)
	{
		sys = XYZSystem(xyzfilepath, size, wanniers);
	}

protected:

	void OutputHeader () const;
};

#endif
