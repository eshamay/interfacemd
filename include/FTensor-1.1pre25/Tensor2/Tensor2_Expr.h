/* Declares a wrapper class for rank 2 Tensor expressions.  I
   specialize it for when I wrap a simple Tensor2 or Tensor2_ptr so
   that it has a reference to the Tensor2(_ptr) and not a copy.
   Otherwise assignment wouldn't work. */

#include "Tensor2_plus_Tensor2.h"
#include "Tensor2_minus_Tensor2.h"
#include "Tensor2_or_Tensor2.h"
#include "Tensor2_times_Tensor2.h"
#include "Tensor2_carat_Tensor2.h"
#include "Tensor2_times_Tensor1.h"
#include "Tensor2_and_Tensor1.h"
#include "Tensor2_times_generic.h"
#include "Tensor2_divide_generic.h"
#include "Tensor2_transform.h"
#include "minus_Tensor2.h"
#include "conj_Tensor2.h"

template<class A, class T, int Dim0, int Dim1, char i, char j>
class Tensor2_Expr
{
  A iter;
public:
  Tensor2_Expr(A &a): iter(a) {}
  T operator()(const int N1, const int N2) const
  {
    return iter(N1,N2);
  }
};

template<class A, class T, int Dim0, int Dim1, char i, char j, Layout layout>
class Tensor2_Expr<Tensor2<A,Dim0,Dim1,layout>,T,Dim0,Dim1,i,j>
{
  Tensor2<A,Dim0,Dim1,layout> &iter;
public:
  Tensor2_Expr(Tensor2<A,Dim0,Dim1,layout> &a): iter(a) {}
  T & operator()(const int N1, const int N2)
  {
    return iter(N1,N2);
  }
  T operator()(const int N1, const int N2) const
  {
    return iter(N1,N2);
  }

  /* Various assignment operators.  I have to explicitly declare the
     second operator= because otherwise the compiler will generate its
     own and not use the template code. */

  template<class B, class U>
  const Tensor2_Expr<Tensor2<A,Dim0,Dim1,layout>,T,Dim0,Dim1,i,j> &
  operator=(const Tensor2_Expr<B,U,Dim0,Dim1,i,j> &result);

  const Tensor2_Expr<Tensor2<A,Dim0,Dim1,layout>,T,Dim0,Dim1,i,j> &
  operator=(const Tensor2_Expr<Tensor2<A,Dim0,Dim1,layout>,T,Dim0,Dim1,i,j> &result);

  template<class B, class U>
  const Tensor2_Expr<Tensor2<A,Dim0,Dim1,layout>,T,Dim0,Dim1,i,j> &
  operator+=(const Tensor2_Expr<B,U,Dim0,Dim1,i,j> &result);

  template<class B, class U>
  const Tensor2_Expr<Tensor2<A,Dim0,Dim1,layout>,T,Dim0,Dim1,i,j> &
  operator-=(const Tensor2_Expr<B,U,Dim0,Dim1,i,j> &result);

  /* This is for when the indices are switched (i,j) -> (j,i). */

  template<class B, class U>
  const Tensor2_Expr<Tensor2<A,Dim0,Dim1,layout>,T,Dim0,Dim1,i,j> &
  operator=(const Tensor2_Expr<B,U,Dim0,Dim1,j,i> &result);

  template<class B, class U>
  const Tensor2_Expr<Tensor2<A,Dim0,Dim1,layout>,T,Dim0,Dim1,i,j> &
  operator+=(const Tensor2_Expr<B,U,Dim0,Dim1,j,i> &result);

  template<class B, class U>
  const Tensor2_Expr<Tensor2<A,Dim0,Dim1,layout>,T,Dim0,Dim1,i,j> &
  operator-=(const Tensor2_Expr<B,U,Dim0,Dim1,j,i> &result);

  /* This is for int's, double's, etc. */

  template <class B>
  const Tensor2_Expr<Tensor2<A,Dim0,Dim1,layout>,T,Dim0,Dim1,i,j> &
  operator=(const B &d);

  template <class B>
  const Tensor2_Expr<Tensor2<A,Dim0,Dim1,layout>,T,Dim0,Dim1,i,j> &
  operator+=(const B &d);

  template <class B>
  const Tensor2_Expr<Tensor2<A,Dim0,Dim1,layout>,T,Dim0,Dim1,i,j> &
  operator-=(const B &d);

  template <class B>
  const Tensor2_Expr<Tensor2<A,Dim0,Dim1,layout>,T,Dim0,Dim1,i,j> &
  operator*=(const B &d);

  template <class B>
  const Tensor2_Expr<Tensor2<A,Dim0,Dim1,layout>,T,Dim0,Dim1,i,j> &
  operator/=(const B &d);
};

/* Specialized for Tensor3_dg_number_rhs_0 (Tensor3_dg with the
   first or second index explicitly given). */

template<class A, class T, int Dim0, int Dim1, char i, char j, int N>
class Tensor2_Expr<Tensor3_dg_number_rhs_0<A,T,N>,T,Dim0,Dim1,i,j>
{
  A &iter;
public:
  Tensor2_Expr(A &a): iter(a) {}
  T & operator()(const int N1, const int N2)
  {
    return iter(N,N1,N2);
  }
  T operator()(const int N1, const int N2) const
  {
    return iter(N,N1,N2);
  }

  /* Various assignment operators.  I have to explicitly declare the
     second operator= because otherwise the compiler will generate its
     own and not use the template code. */

  template<class B, class U>
  const Tensor2_Expr<Tensor3_dg_number_rhs_0<A,T,N>,T,Dim0,Dim1,i,j> &
  operator=(const Tensor2_Expr<B,U,Dim0,Dim1,i,j> &result);

  const Tensor2_Expr<Tensor3_dg_number_rhs_0<A,T,N>,T,Dim0,Dim1,i,j> &
  operator=(const Tensor2_Expr<Tensor3_dg_number_rhs_0
	    <A,T,N>,T,Dim0,Dim1,i,j> &result);

};

#include "Tensor2_Expr_equals.h"
