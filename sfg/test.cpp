#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>

#include "../matrixr.h"

using namespace boost::numeric::ublas;
using std::cout;
using std::endl;


void Print (const symmetric_matrix<double>& t) {
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

  symmetric_matrix<double> m(5,5);
  m(1,4) = 5.0;
  m(2,1) = 2.0;
  m.clear();

  Print(m);
  return 0;
}
