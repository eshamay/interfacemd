#pragma once
#ifndef AMBER_MORITA2002_H_
#define AMBER_MORITA2002_H_

#include "morita2002.h"

namespace morita {
  class CP2KMorita2002Analysis : public Morita2002Analysis<XYZSystem> {
	public:

	  CP2KMorita2002Analysis (std::string fn) : 
		Morita2002Analysis<XYZSystem>(
			std::string("Analysis of the waters in an MD system using the technique of Morita/Hynes (2002) with a system produced by the CP2K package"),
			fn) { }

	protected:
	  void SelectAnalysisWaters ();
	  void SetAnalysisWaterDipoleMoments ();

  };	// class CP2K sfg analyzer

} // namespace morita

#endif

