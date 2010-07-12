#include "mdmath.h"


/*
namespace scalapack {



  PDGEMM::PDGEMM (tensor::tensor_t& A, tensor::tensor_t& B, tensor::tensor_t& C, boost::mpi::communicator& mpi_comm)
	:
	  world(mpi_comm),
	  nrow(2), ncol(world.size()/nrow), nblock(64),
	  ma(A.size1()), mb(B.size1()), mc(C.size1()),
	  na(A.size2()), nb(B.size2()), nc(C.size2()),
	  transa('N'), transb('N'),
	  alpha(1.0), beta(1.0),
	  scope('A'), top(' ')
  {

	// initialize the processor grid
	blacs::sl_init_ (&ctxt, &nrow, &ncol);
	blacs::blacs_gridinfo_ (&ctxt, &nrow, &ncol, &lnrow, &lncol);

	// determine local matrix dimensions
	lma = numroc_ (&ma, &nblock, &lnrow, &izero, &nrow);
	lna = numroc_ (&na, &nblock, &lncol, &izero, &ncol);	
	tensor::tensor_t lA (lma, lna);

	lmb = numroc_ (&mb, &nblock, &lnrow, &izero, &nrow);
	lnb = numroc_ (&nb, &nblock, &lncol, &izero, &nrow);	
	tensor::tensor_t lB (lmb, lnb);

	lmc = numroc_ (&mc, &nblock, &lnrow, &izero, &nrow);
	lnc = numroc_ (&nc, &nblock, &lncol, &izero, &nrow);	
	tensor::tensor_t lC (lmc, lnc);

	// distribute the global matrices to the processes
	boost::mpi::broadcast(world,&A(0,0),ma*mb,0);
	boost::mpi::broadcast(world,&B(0,0),na*nb,0);

	// set up the matrix descriptors for the scalapack routine
	descinit_ (desca, &ma, &na, &nblock, &nblock, &izero, &izero, &ctxt, &lma, &ierr);
	descinit_ (descb, &mb, &nb, &nblock, &nblock, &izero, &izero, &ctxt, &lmb, &ierr);
	descinit_ (descc, &mc, &nc, &nblock, &nblock, &izero, &izero, &ctxt, &lmc, &ierr);

	// 2d block-cyclic data distribution to the processes
	double val;
	for (int i = 1; i <= ma; i++) {
	  for (int j = 1; j <= na; j++) {
		val = A(i-1,j-1);
		pdelset_ (&lA(0,0), &i, &j, desca, &val);
	  }
	}
	for (int i = 1; i <= mb; i++) {
	  for (int j = 1; j <= nb; j++) {
		val = B(i-1,j-1);
		pdelset_ (&lB(0,0), &i, &j, descb, &val);
	  }
	}

	// perform the local solution of the multiplication
	pdgemm_ (&transa, &transb, &ma, &nb, &na, &alpha, &lA(0,0), &ione, &ione, desca, &lB(0,0), &ione, &ione, descb, &beta, &lC(0,0), &ione, &ione, descc);

	// gather the local solutions to the global one
	for (int i = 1; i <= mc; i++) {
	  for (int j = 1; j <= nc; j++) {
		pdelget_ (&scope, &top, &C(i-1,j-1), &lC(0,0), &i, &j, descc);
	  }
	}

	blacs::blacs_gridexit_ (&ctxt);

  } // pdgemm c-tor

}	// namespace scalapack
*/


