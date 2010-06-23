#include "tensor.h"

namespace tensor {
  using namespace boost::numeric::ublas;

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
	  Tensor<T> m;
	  m = boost::numeric::ublas::trans(*this);
	  return m;
	}


  template <typename T>
	double Tensor<T>::Trace () const {
	  return (*this)(0,0) + (*this)(1,1) + (*this)(2,2);
	} //trace


  template <typename T>
	double Tensor<T>::Determinant () const
	{
	  double det = 1.0;

	  matrix_t mLu (*this);
	  permutation_matrix<std::size_t> pivots(this->size1() );

	  int is_singular = lu_factorize(mLu, pivots);

	  if (!is_singular)
	  {
		for (std::size_t i=0; i < pivots.size(); ++i)
		{
		  if (pivots(i) != i)
			det *= -1.0;

		  det *= mLu(i,i);
		}
	  }
	  else
		det = 0.0;

	  return det;

	} // determinant


} // namespace tensor
