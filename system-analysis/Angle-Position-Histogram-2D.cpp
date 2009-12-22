#include "Angle-Position-Histogram-2D.h"

int main () {
  WaterSystemParams wsp ("2D-Angle-Position-Histogram.Normal.dat");

  Angle_Position_Analyzer apa (wsp);

  apa.SystemAnalysis();
  return 0;
}


