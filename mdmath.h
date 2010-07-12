#ifndef MDMATH_H_
#define MDMATH_H_

#include "tensor.h"
#include <mkl_blas.h>
#include <mkl_lapack.h>

/*
namespace blacs {

  extern "C" {
	void sl_init_			(int* ictxt, int* nprow, int* npcol);
	void blacs_gridinfo_	(int* context, int *nprow, int *npcol, int *myrow, int *mycol);
	void blacs_pinfo_		(int* id, int* nprocs);
	void blacs_get_			(int* context, int* request, int* value);
	int  blacs_gridinit_	(int* context, char * order, int* np_row, int* np_col);
	void blacs_gridexit_	(int* context);
	void blacs_exit_		(int* error_code);


	   void   Cblacs_pinfo( int* mypnum, int* nprocs);
	   void   Cblacs_get( int context, int request, int* value);
	   int    Cblacs_gridinit( int* context, char * order, int np_row, int np_col);
	   void   Cblacs_gridinfo( int context, int*  np_row, int* np_col, int*  my_row, int*  my_col);
	   void   Cblacs_gridexit( int context);
	   void   Cblacs_exit( int error_code);

	   double pdlamch_( int *ictxt , char *cmach);
	   double pdlange_( char *norm, int *m, int *n, double *A, int *ia, int *ja, int *desca, double *work);
  }	// extern

} // blacs
*/

/*
namespace lapack {

  extern "C" {

	void dgesv_ (int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);

  } // extern

} // namespace lapack
*/



/*
namespace scalapack {

  extern "C" {
	int    numroc_( int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);

	void   descinit_( int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc, int *ictxt, int *lld, int *info); 

	void pdgemm_( char *TRANSA, char *TRANSB, int * M, int * N, int * K, double * ALPHA, double * A, int * IA, int * JA, int * DESCA, double * B, int * IB, int * JB, int * DESCB, double * BETA, double * C, int * IC, int * JC, int * DESCC );

	void pdelset_ (double * array, int * row, int * col, int * desca, double * value);
	void pdelget_ (char * scope, char * top, double * alpha, double * a, int * ia, int * ja, int * desca);
  } // extern



} // scalapack
*/

#endif
