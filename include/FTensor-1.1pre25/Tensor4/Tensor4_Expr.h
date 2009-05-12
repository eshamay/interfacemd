/* Declare a wrapper class for generic rank 4 Tensor expressions. */

#include "Tensor4_plus_Tensor4.h"
#include "Tensor4_times_Tensor2_symmetric.h"
#include "Tensor4_times_Tensor2.h"

template<class A, class T, int Dim0, int Dim1, int Dim2, int Dim3,
  char i, char j, char k, char l>
class Tensor4_Expr
{
  A iter;
public:
  Tensor4_Expr(A &a): iter(a) {}
  T operator()(const int N1, const int N2, const int N3,
	       const int N4) const
  {
    return iter(N1,N2,N3,N4);
  }
};
