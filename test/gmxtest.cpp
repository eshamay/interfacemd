#include "gmxtest.h"

void GMXAnalyzer::Setup () { 
  FindWaters("SM2");
  return;
}

void GMXAnalyzer::Analysis () {

  return;
}

// nothing at the moment
void GMXAnalyzer::PostAnalysis () { return; }

void GMXAnalyzer::DataOutput (const unsigned int timestep) {
  //if (!(timestep % (output_freq * 10)))
  //binner.Output(output, timestep);
  return;
}


int main () {

  libconfig::Config cfg;
  cfg.readFile("system.cfg");

  WaterSystemParams wsp (cfg);

  GMXAnalyzer an (wsp);

  an.SystemAnalysis();


  return 0;

}
