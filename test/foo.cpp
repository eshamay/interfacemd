#include <boost/timer.hpp>
#include "../tensor.h"

#include <Eigen/Core>
  USING_PART_OF_NAMESPACE_EIGEN

using namespace std;

int main () {

  Matrix<double,12,2> m;
  Matrix<double,3,2> v;
  v << 1,2,3,4,5,6;

  cout << m << endl;
  cout << endl;
  cout << v << endl;

  for (int i = 0; i < 4; i++) {
	m.block(3*i,0,3,2) = v;
  }

  cout << endl;
  cout << m << endl;

  return 0;
}
