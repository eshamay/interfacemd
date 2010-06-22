#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>
#include "../vecr.h"

//#include "../matrixr.h"

using namespace boost::numeric::ublas;
using std::cout;
using std::endl;


void Print (const matrix<double>& t) {
  for (unsigned i = 0; i < t.size1(); i++) {
	for (unsigned j = 0; j < t.size2(); j++) {
	  printf ("% 8.3f", t(i,j));
	}
	printf ("\n");
  }
  printf ("\n");
}

int main () 
{

  matrix<double> m (4,3);

  VecR u (2.0, 5.0, 3.0);
  VecR v (u);
  VecR w;
  double p[3] = {4.0, 3.0, 7.0};
  VecR z(p);

  cout << u << endl;
  u.Zero();
  cout << u << endl;

  return 0;
}
