#include "angle-histogram.h"

AngleBinner::histo_t AngleBinner::histo = AngleBinner::histo_t (-1.0, 1.0, 0.01);

int main () {

  WaterSystemParams wsp ("angle-histogram.dat", 75000);

  AngleHistogram ah (wsp);

  ah.SystemAnalysis();

  return 0;
}
