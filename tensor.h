#ifndef __TENSOR_H
#define __TENSOR_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/lu.hpp>

namespace tensor {

  using namespace boost::numeric::ublas;

  typedef matrix<double> matrix_t;
  typedef identity_matrix<double> id_matrix_t;

  template <typename T>
	class Tensor : public T {

	  public:
		Tensor (const int M, const int N) : T(M, N) { }
		Tensor (const matrix_t& t) : T(t) { }
		virtual ~Tensor() { }

		virtual void Print() const;

		void 	Zero ();
		double 	Determinant () const;
		double	Trace () const;
		Tensor<T>	Transpose () const;


	}; // tensor

  typedef Tensor<matrix_t>	tensor_t;
  typedef Tensor<symmetric_matrix<double, upper> > sym_tensor_t;



  // A matrix that forms an Mx3 tensor where each 3x3 submatrix is the identity matrix
  class ColumnIdentityMatrix : public tensor_t {

	public:
	  ColumnIdentityMatrix (const int M) : tensor_t(M,3)
	{
	  id_matrix_t id (3,3);
	  for (unsigned i = 0; i < size1()/3; i++)
		project(*this, slice(3*i,1,3), slice(0,1,3)) = id;
	}
  }; // Column Identity Matrix



  class SymmetricMatrix : public sym_tensor_t {
	public:
	  SymmetricMatrix () : sym_tensor_t(0,0) { }
  }; // symmetric matrix


} // namespace tensor
#endif
