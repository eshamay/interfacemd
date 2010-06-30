#include "test.h"

using std::cout;
using std::endl;

  int main (int argc, char **argv) 
  {

	// set up the mpi stack
	boost::mpi::environment env(argc, argv);
	boost::mpi::communicator world;

	int n = 900;
	tensor::tensor_t A(n,n);
	tensor::tensor_t B(n,n);
	tensor::tensor_t C(n,n);
	tensor::tensor_t D(n,n);
	tensor::tensor_t E(n,n);


	if (world.rank() == 0) {
	  // load the data into the matrices
	  for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
		  A(i,j) = i+j;
		  B(i,j) = 3*i+j;
		}
	  }
	}

	scalapack::PDGEMM p (A,B,C, world);


	/*
	if (!world.rank()) {
	  D.assign(prod(A,B));
	  E.assign(D-C);
	  bool set = true;
	  for (int i = 0; i < n && set; i++) {
		for (int j = 0; j < n && set; j++) {
		  if (E(i,j) != 0.0) printf ("doh!\n");
		  set = false;
		}
	  }
	}
	*/

	/*

	int mpi_size = world.size();
	int mpi_rank = world.rank();

	// set up the processor grid for blacs
	int nr = 2, nc = mpi_size/nr;	// blacs grid
	int ctxt;	// blacs context

	blacs::sl_init_ (&ctxt, &nr, &nc);

	int lnr, lnc;	// local grid assignment
	blacs::blacs_gridinfo_ (&ctxt, &nr, &nc, &lnr, &lnc);



	// Figure out A*B = C. Create A and B, find C
	// data to operate on
	int n = 128;
	tensor::tensor_t A(n,n);
	tensor::tensor_t B(n,n);
	tensor::tensor_t C(n,n);

	if (mpi_rank == 0) {
	  // load the data into the matrices
	  for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
		  A(i,j) = i+j;
		  B(i,j) = n*i+j;
		}
	  }
	  //A.Print();
	}


	int nb = 4;	// matrix block size for cyclic distribution
	int izero = 0, ione = 1;
	// determine the local array sizes for the computations
	int lna = blacs::numroc_ (&n, &nb, &lnr, &izero, &nr);	// local matrix sizes
	int lnb = blacs::numroc_ (&n, &nb, &lnc, &izero, &nc);


	int desca[9];
	int descb[9];
	int descc[9];
	int ierr;

	blacs::descinit_ (desca, &n, &n, &nb, &nb, &izero, &izero, &ctxt, &lna, &ierr);
	blacs::descinit_ (descb, &n, &n, &nb, &nb, &izero, &izero, &ctxt, &lna, &ierr);
	blacs::descinit_ (descc, &n, &n, &nb, &nb, &izero, &izero, &ctxt, &lna, &ierr);

	//tensor::tensor_t lA(lna,lnb);
	tensor::tensor_t lA(lna,lnb);
	tensor::tensor_t lB(lna,lnb);
	tensor::tensor_t lC(lna,lnb);


	// distribute the global matrices to each process
	boost::mpi::broadcast(world,&A(0,0),n*n,0);
	boost::mpi::broadcast(world,&B(0,0),n*n,0);

	double val;
	for (int i = 1; i <= n; i++) {
	  for (int j = 1; j <= n; j++) {
		val = A(i-1,j-1);
		blacs::pdelset_ (&lA(0,0), &i, &j, desca, &val);
		val = B(i-1,j-1);
		blacs::pdelset_ (&lB(0,0), &i, &j, descb, &val);
	  }
	}

	char transa = 'N';
	char transb = 'N';
	double alpha = 1.0, beta = 1.0;
	// matrix multiplication
	blacs::pdgemm_ (&transa, &transb, &n, &n, &n, &alpha, &lA(0,0), &ione, &ione, desca, &lB(0,0), &ione, &ione, descb, &beta, &lC(0,0), &ione, &ione, descc);

	char scope = 'A', top = ' ';
	for (int i = 1; i <= n; i++) {
	  for (int j = 1; j <= n; j++) {
		blacs::pdelget_ (&scope, &top, &C(i-1,j-1), &lC(0,0), &i, &j, descc);
	  }
	}


	tensor::tensor_t D(n,n);
	tensor::tensor_t E(n,n);
	D.assign(prod(A,B));
	E.assign(C-D);
	for (int i = 0; i < n; i++) {
	  for (int j = 0; j < n; j++) {
		if (E(i,j)) {
		  printf ("wrong!\n");
		  break;
		}
	  }
	}

	blacs::blacs_gridexit_ (&ctxt);
	*/

	return 0;
  }
