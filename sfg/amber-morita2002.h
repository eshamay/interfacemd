#pragma once
#ifndef AMBER_MORITA2002_H_
#define AMBER_MORITA2002_H_

#include "morita2002.h"

class AmberSFGAnalyzer : public SFGAnalyzer<AmberSystem> {

public:
  AmberSFGAnalyzer (WaterSystemParams& wsp);

protected:
  void SelectAnalysisWaters ();
  void SetAnalysisWaterDipoleMoments ();

};	// class Amber sfg analyzer

#endif
