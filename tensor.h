#ifndef __TENSOR_H
#define __TENSOR_H

#define EIGEN_MATRIXBASE_PLUGIN "/home/eshamay/md/src/EigenMatrixAddon.h"
#include <Eigen/Core>
#include <list>

namespace tensor {

  USING_PART_OF_NAMESPACE_EIGEN

  // Turn every sub-block (of given blocksize = bs) within the matrix into an identity matrix
  template <typename Derived>
  void BlockIdentity (MatrixBase<Derived>& m, int bs) {

	MatrixXd block_m (bs,bs);
	block_m.setIdentity();

	for (unsigned i = 0; i < m.rows()/bs; i++) {
	  for (unsigned j = 0; j < m.cols()/bs; j++) {
		m.block(i*bs,j*bs,bs,bs) = block_m;
		//project (m, slice(bs*i,1,bs), slice(bs*j,1,bs)) = id;
	  } 
	}
  }	// block identity

} // namespace tensor
#endif
