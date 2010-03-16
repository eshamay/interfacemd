#include "Angle-Position-Histogram-2D.h"

int main () {
  libconfig::Config cfg;
  cfg.readFile("system.cfg");

  std::string filename = cfg.lookup("analysis.angle-position-histogram.filename");
  libconfig::Setting &analysis = cfg.lookup("analysis");
  analysis.add("filename", libconfig::Setting::TypeString) = filename;

  //WaterSystemParams wsp ("2D-Angle-Position-Histogram.Normal.dat", 100000);
  //WaterSystemParams wsp (cfg, "2D-Angle-Position-Histogram.Bisector.dat", 100000);
  WaterSystemParams wsp(cfg);

  Angle_Position_Analyzer apa (wsp);
  apa.SystemAnalysis();



  return 0;
}


