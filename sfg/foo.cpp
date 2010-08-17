#include "foo.h"

int main () {

  // for now, f is 'h', the 3Nx3 block identity tensor
  MatR A;
  A << 1,2,3,8,7,1,2,2,1;
  A.Print();
  //std::vector<double> i = A.ColumnMajorMatrixToVector();
  //std::copy(&A(0,0), &A(A.rows(),A.cols()-1), std::ostream_iterator<double>(std::cout, "\n"));
  std::ofstream fout ("foo.dat");
  std::copy(&A(0,0), &A(A.rows(),A.cols()-1), std::ostream_iterator<double>(fout, "\n"));
  
  for (unsigned int i = 0; i < 3; i++) {
	for (unsigned int j = 0; j < 3; j++) {
	  printf ("% 8.3f", A(i,j));
	}
  }
  printf ("\n");
  fout.close();


}
