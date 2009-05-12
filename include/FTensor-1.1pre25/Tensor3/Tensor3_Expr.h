/* Declare a wrapper class for generic rank 3 Tensor expressions.
   There isn't a Tensor3 class yet.  I only use Tensor3_Expr as an
   intermediate expression which immediately get contracted with
   something to make a Tensor2 or Tensor1. */

#include "Tensor3_times_generic.h"
#include "Tensor3_times_Tensor1.h"
#include "Tensor3_times_Tensor2_symmetric.h"
#include "Tensor3_times_Tensor2.h"
#include "Tensor3_times_Tensor3.h"
#include "Tensor3_times_Tensor3_dg.h"
#include "Tensor3_plus_Tensor3.h"
#include "Tensor3_or_Tensor3.h"
#include "Tensor3_minus_Tensor3_dg.h"

template<class A, class T, int Dim0, int Dim1, int Dim2,
  char i, char j, char k>
class Tensor3_Expr
{
  A iter;
public:
  Tensor3_Expr(A &a): iter(a) {}
  T operator()(const int N1, const int N2, const int N3) const
  {
    return iter(N1,N2,N3);
  }
};
