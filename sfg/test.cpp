#include "test.h"

using std::cout;
using std::endl;

int main (int argc, char **argv) 
{

  tensor::tensor_t A (3,3);
  tensor::tensor_t B (3,3);
  tensor::tensor_t C (3,3);

  for (int i = 0; i < 3; i++) {
	for (int j = 0; j < 3; j++) {
	  A(i,j) = 3*i+j;
	  B(i,j) = i+j;
	}
  }

  A.Print();
  B.Print();

  C.assign(prod(A,B));
  C.Print();

  return 0;
}
