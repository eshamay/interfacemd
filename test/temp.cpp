#include <iostream>
#include "../analysis.h"

using namespace std;

typedef pair<int,int> pair_t;

int main () {

  WaterSystemParams wsp ("temp.out");
  Analyzer an(wsp);

  return 0;
}
