#include "../matrixr.h"

using std::cout;
using std::endl;

int main () 
{

  double dm[9] = {3,2,1,9,7,6,5,8,4};
  MatR m(dm);
  double dn[9] = {8,4,5,6,3,1,9,2,7};
  MatR n(dn);

  m.Print();
  n.Print();

  (m + 3.0).Print();


  return 0;
}
