#ifndef __TENSOR_H
#define __TENSOR_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>

namespace tensor {

  using namespace boost::numeric::ublas;

  typedef matrix<double> matrix_t;
  typedef identity_matrix<double> id_matrix_t;

  template <typename T>
	class Tensor : public T {

	  public:
		Tensor (const int M, const int N) : T(M, N) { }
		Tensor (const Tensor<T>& t) : T(t) { }
		virtual ~Tensor() { }


		void 	Zero ();
		int determinant_sign(const permutation_matrix<std ::size_t>& pm) const;
		double 	Determinant () const;
		void 	Inverse (T& inv) const;
		double	Trace () const;
		Tensor<T>	Transpose () const;

		void Print() const;


	}; // tensor
  typedef Tensor<matrix_t>	tensor_t;
  typedef Tensor<symmetric_matrix<double, upper> > sym_tensor_t;


  template <typename T>
	void Tensor<T>::Print() const
	{ 
	  for (unsigned i = 0; i < this->size1(); i++) {
		for (unsigned j = 0; j < this->size2(); j++) {
		  printf ("% 8.3f", (*this)(i,j));
		}
		printf ("\n");
	  }
	  printf ("\n");
	}	// Print

  template <typename T>
	void Tensor<T>::Zero() {
	  for (unsigned i = 0; i < this->size1(); i++) {
		for (unsigned j = 0; j < this->size2(); j++) {
		  (*this)(i,j) = T(0);
		} 
	  }
	}

  template <typename T>
	Tensor<T> Tensor<T>::Transpose () const {
	  //Tensor<T> m (this->size1(), this->size2());
	  Tensor<T> m (*this);
	  m.assign (trans(*this));
	  return m;
	}


  template <typename T>
	double Tensor<T>::Trace () const {
	  return (*this)(0,0) + (*this)(1,1) + (*this)(2,2);
	} //trace




  /** General matrix inversion routine.
   * It uses lu_factorize and lu_substitute in uBLAS to invert a matrix
   */
  template<class T>
	void Tensor<T>::Inverse (T& inv) const {
	  // create a working copy of the input
	  matrix_t mLu(*this);

	  // perform LU-factorization
	  lu_factorize(mLu);

	  // create identity matrix of "inverse"
	  inv.assign(id_matrix_t(mLu.size1()));

	  // backsubstitute to get the inverse
	  lu_substitute<matrix_t const, T >(mLu, inv);

	} // lu_inv



  template <typename T>
	int Tensor<T>::determinant_sign(const permutation_matrix<std ::size_t>& pm) const
	{
	  int pm_sign=1;
	  std::size_t size = pm.size();
	  for (std::size_t i = 0; i < size; ++i)
		if (i != pm(i))
		  pm_sign *= -1; // swap_rows would swap a pair of rows here, so we change sign
	  return pm_sign;
	}



  template <typename T>
	double Tensor<T>::Determinant() const 
	{

	  matrix_t dm (*this);
	  permutation_matrix<std::size_t> pm (dm.size1());
	  double det = 1.0;
	  if( lu_factorize(dm,pm) ) {
		det = 0.0;
	  } else {
		for (unsigned i = 0; i < dm.size1(); i++) 
		  det *= dm(i,i); // multiply by elements on diagonal
		det = det * determinant_sign( pm );
	  }
	  return det;

	} // Determinant




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
