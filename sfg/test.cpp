#include "../matrixr.h"

using std::cout;
using std::endl;

int main () 
{

  double dm[9] = {2,-1,3,-1,-1,3,1,-2,-1};
  MatR m(dm);
  m.Print();
  //double dn[9] = {8,4,5,6,3,1,9,2,7};
  //MatR n(dn);
  //n.Print();

  VecR v (-3,-6,-2);

  MatR inv;
  tensor::Tensor<tensor::matrix_t>::Inverse(m,inv);
  (inv*v).Print();

  m.Zero();
  m.Print();
  cout << m.size1() << " " << m.size2() << endl;

  return 0;
}
