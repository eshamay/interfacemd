#pragma once
#ifndef AMBER_MORITA2002_H_
#define AMBER_MORITA2002_H_

#include "morita2002.h"

namespace morita {
  class CP2KSFGAnalyzer : public SFGAnalyzer<XYZSystem> {
	public:

	  CP2KSFGAnalyzer (WaterSystemParams& wsp);

	protected:
	  void SelectAnalysisWaters ();
	  void SetAnalysisWaterDipoleMoments ();

  };	// class CP2K sfg analyzer

} // namespace morita

#endif

