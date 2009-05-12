/* Declares a wrapper class for rank 4 Tensor expressions symmetric on
   the first two and last two indices. */

#include "Tensor4_ddg_times_Tensor2_symmetric.h"
#include "Tensor4_ddg_carat_Tensor2_symmetric.h"
#include "Tensor4_ddg_and_Tensor2_symmetric.h"
#include "Tensor4_ddg_mod_Tensor2_symmetric.h"
#include "Tensor4_ddg_times_Tensor2.h"
#include "Tensor4_ddg_times_Tensor1.h"
#include "Tensor4_ddg_times_generic.h"
#include "Tensor4_ddg_times_Tensor4_ddg.h"
#include "Tensor4_ddg_plus_Tensor4_ddg.h"
#include "Tensor4_ddg_minus_Tensor4_ddg.h"
#include "Tensor4_ddg_or_Tensor4_ddg.h"
#include "Tensor4_ddg_and_Tensor4_ddg.h"
#include "Tensor4_ddg_carat_Tensor4_ddg.h"
#include "minus_Tensor4_ddg.h"
//  #include "Tensor4_ddg_mod_Tensor4_ddg.h"

template<class A, class T, int Dim01, int Dim23,
  char i, char j, char k, char l>
class Tensor4_ddg_Expr
{
  A iter;
public:
  Tensor4_ddg_Expr(A &a): iter(a) {}
  T operator()(const int N1, const int N2, const int N3, const int N4)
    const
  {
    return iter(N1,N2,N3,N4);
  }
};

template<class A, class T, int Tensor_Dim01, int Tensor_Dim23,
  int Dim01, int Dim23, char i, char j, char k, char l>
class Tensor4_ddg_Expr<Tensor4_ddg<A,Tensor_Dim01,Tensor_Dim23>,T,Dim01,Dim23,i,j,k,l>
{
  Tensor4_ddg<A,Tensor_Dim01,Tensor_Dim23> &iter;
public:
  Tensor4_ddg_Expr(Tensor4_ddg<A,Tensor_Dim01,Tensor_Dim23> &a): iter(a) {}
  T operator()(const int N1, const int N2, const int N3, const int N4)
    const
  {
    return iter(N1,N2,N3,N4);
  }

  /* Various assignment operators.  I have to explicitly declare the
     second operator= because otherwise the compiler will generate its
     own and not use the template code. */

  template<class B, class U>
  const Tensor4_ddg_Expr<Tensor4_ddg<A,Tensor_Dim01,Tensor_Dim23>,T,Dim01,Dim23,i,j,k,l> &
  operator=(const Tensor4_ddg_Expr<B,U,Dim01,Dim23,i,j,k,l> &result);

  const Tensor4_ddg_Expr<Tensor4_ddg<A,Tensor_Dim01,Tensor_Dim23>,T,Dim01,Dim23,i,j,k,l> &
  operator=(const Tensor4_ddg_Expr<Tensor4_ddg<A,Tensor_Dim01,Tensor_Dim23>,T,Dim01,Dim23,i,j,k,l> &result);
};

#include "Tensor4_ddg_Expr_equals.h"
