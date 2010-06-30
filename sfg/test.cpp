#include "test.h"

using std::cout;
using std::endl;

  int main (int argc, char **argv) 
  {

	// set up the mpi stack
	boost::mpi::environment env(argc, argv);
	boost::mpi::communicator world;

	int n = 512;
	tensor::tensor_t A(n,n);
	tensor::tensor_t B(n,n);
	tensor::tensor_t C(n,n);
	tensor::tensor_t D(n,n);


	if (world.rank() == 0) {
	  // load the data into the matrices
	  for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
		  A(i,j) = i+j;
		  B(i,j) = 3*i+j;
		}
	  }
	}

	char transa = 'N';
	char transb = 'N';
	double scale = 1.0;

	for (int i = 0; i < 128; i++) {
	  blas::dgemm_ (&transa, &transb, &n, &n, &n, &scale, &A(0,0), &n, &B(0,0), &n, &scale, &C(0,0), &n);
	  //scalapack::PDGEMM p (A,B,C, world);
	}



	return 0;
  }
