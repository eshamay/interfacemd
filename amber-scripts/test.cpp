#include "quartz-utils.h"

using namespace std;

int main () {

  int base = 26;
  cout << itoa_base (0,base) << endl;
  cout << itoa_base (1,base) << endl;
  cout << itoa_base (25,base) << endl;
  cout << itoa_base (192,base) << endl;

  return 0;
}
