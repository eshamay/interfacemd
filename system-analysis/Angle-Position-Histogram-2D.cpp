#include "Angle-Position-Histogram-2D.h"

int main () {
  WaterSystemParams wsp ("2D-Angle-Position-Histogram.Normal.dat", 75000);
  //WaterSystemParams wsp ("2D-Angle-Position-Histogram.Bisector.dat", 75000);

  Angle_Position_Analyzer apa (wsp);

  apa.SystemAnalysis();
  return 0;
}


