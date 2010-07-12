#include "test.h"

using std::cout;
using std::endl;

int main (int argc, char **argv) 
{

  int N = 1024;
  tensor::tensor_t A (N,N);
  tensor::tensor_t B (N,N);
  tensor::tensor_t C (N,N);
  tensor::tensor_t D (N,N);

  for (int i = 0; i < N; i++) {
	for (int j = 0; j < N; j++) {
	  A(i,j) = N*i+j;
	  B(i,j) = i+j;
	}
  }

  //A.Print();
  //B.Print();

  //D.assign(prod(A,B));
  //C.Print();

  C.clear();
  char trans = 'N';
  double scale = 1.0;

  dgemm (&trans, &trans, &N, &N, &N, &scale, &A(0,0), &N, &B(0,0), &N, &scale, &C(0,0), &N);
  //C.Print();


  /*
  for (int i = 0; i < N; i++) {
	for (int j = 0; j < N; j++) {
	  if (C(i,j) != D(i,j)) {
		printf ("WRONG!\n");
	  } 
	}
  }
  */



  return 0;
}
