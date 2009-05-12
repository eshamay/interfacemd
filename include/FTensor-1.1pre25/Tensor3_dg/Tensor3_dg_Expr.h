/* Declares a wrapper class for symmetric rank 3 Tensor expressions
   that are symmetric on the first two indices.  It is also used for
   Tensor3_christof's but with some games played with the
   indices. */

#include "Tensor3_dg_plus_Tensor3_dg.h"
#include "Tensor3_dg_minus_Tensor3_dg.h"
#include "Tensor3_dg_and_Tensor3_dg.h"
#include "Tensor3_dg_or_Tensor3_dg.h"
#include "Tensor3_dg_times_Tensor3_dg.h"
#include "Tensor3_dg_times_Tensor2.h"
#include "Tensor3_dg_times_Tensor2_symmetric.h"
#include "Tensor3_dg_and_Tensor2_symmetric.h"
#include "Tensor3_dg_times_Tensor1.h"
#include "Tensor3_dg_and_Tensor1.h"
#include "Tensor3_dg_times_generic.h"
#include "Tensor3_dg_divide_generic.h"
#include "minus_Tensor3_dg.h"

template<class A, class T, int Dim01, int Dim2, char i, char j, char k>
class Tensor3_dg_Expr
{
  A iter;
public:
  Tensor3_dg_Expr(A &a): iter(a) {}
  T operator()(const int N1, const int N2, const int N3) const
  {
    return iter(N1,N2,N3);
  }
};

template<class A, class T, int Tensor_Dim01, int Tensor_Dim2,
  int Dim01, int Dim2, char i, char j, char k>
class Tensor3_dg_Expr<Tensor3_dg<A,Tensor_Dim01,Tensor_Dim2>,T,Dim01,Dim2,i,j,k>
{
  Tensor3_dg<A,Tensor_Dim01,Tensor_Dim2> &iter;
public:
  Tensor3_dg_Expr(Tensor3_dg<A,Tensor_Dim01,Tensor_Dim2> &a): iter(a) {}
  T & operator()(const int N1, const int N2, const int N3)
  {
    return iter(N1,N2,N3);
  }
  T operator()(const int N1, const int N2, const int N3) const
  {
    return iter(N1,N2,N3);
  }

  /* Various assignment operators.  I have to explicitly declare the
     second operator= because otherwise the compiler will generate its
     own and not use the template code. */

  template<class B, class U> inline
  const Tensor3_dg_Expr<Tensor3_dg<A,Tensor_Dim01,Tensor_Dim2>,T,Dim01,Dim2,i,j,k> &
  operator=(const Tensor3_dg_Expr<B,U,Dim01,Dim2,i,j,k> &result);

  const Tensor3_dg_Expr<Tensor3_dg<A,Tensor_Dim01,Tensor_Dim2>,T,Dim01,Dim2,i,j,k> &
  operator=(const Tensor3_dg_Expr<Tensor3_dg<A,Tensor_Dim01,Tensor_Dim2>,T,Dim01,Dim2,i,j,k> &result);

  template <class U> inline
  const Tensor3_dg_Expr<Tensor3_dg<A,Tensor_Dim01,Tensor_Dim2>,T,Dim01,Dim2,i,j,k> &
  operator=(const U &d);
};

/* I need a version for const and non-const Tensor3_christof,
   otherwise it will use the default and the index order will get all
   messed up. The class A is either T or T*, depending on whether
   Tensor3_christof has an array of T's or T *'s. */

template<class A, class T, int Tensor_Dim0, int Tensor_Dim12,
  int Dim12, int Dim0, char i, char j, char k>
class Tensor3_dg_Expr<const Tensor3_christof<A,Tensor_Dim0,Tensor_Dim12>,
  T,Dim12,Dim0,i,j,k>
{
  const Tensor3_christof<A,Tensor_Dim0,Tensor_Dim12> &iter;
public:
  Tensor3_dg_Expr(const Tensor3_christof<A,Tensor_Dim0,Tensor_Dim12> &a): iter(a) {}

  /* Need to switch the index order because a christof tensor is just
     a dg tensor with the indices switched. */

  T operator()(const int N1, const int N2, const int N3) const
  {
    return iter(N3,N1,N2);
  }
};

template<class A, class T, int Tensor_Dim0, int Tensor_Dim12,
  int Dim12, int Dim0, char i, char j, char k>
class Tensor3_dg_Expr<Tensor3_christof<A,Tensor_Dim0,Tensor_Dim12>,T,Dim12,Dim0,i,j,k>
{
  Tensor3_christof<A,Tensor_Dim0,Tensor_Dim12> &iter;
public:
  Tensor3_dg_Expr(Tensor3_christof<A,Tensor_Dim0,Tensor_Dim12> &a): iter(a) {}

  /* Need to switch the index order because a christof tensor is just
     a dg tensor with the indices switched.  I have to explicitly
     declare the second operator= because otherwise the compiler will
     generate its own and not use the template code. */

  T operator()(const int N1, const int N2, const int N3) const
  {
    return iter(N3,N1,N2);
  }
  template<class B, class U> inline
  const Tensor3_dg_Expr<Tensor3_christof<A,Tensor_Dim0,Tensor_Dim12>,T,Dim12,Dim0,i,j,k> &
  operator=(const Tensor3_dg_Expr<B,U,Dim12,Dim0,i,j,k> &result);

  inline
  const Tensor3_dg_Expr<Tensor3_christof<A,Tensor_Dim0,Tensor_Dim12>,T,Dim12,Dim0,i,j,k> &
  operator=(const Tensor3_dg_Expr<Tensor3_christof<A,Tensor_Dim0,Tensor_Dim12>,
	    T,Dim12,Dim0,i,j,k> &result);

  template<class U> inline
  const Tensor3_dg_Expr<Tensor3_christof<A,Tensor_Dim0,Tensor_Dim12>,T,Dim12,Dim0,i,j,k> &
  operator=(const U &d);
};

/* Specialized for Tensor4_ddg_number_rhs_0 (Tensor4_ddg with the
   first index explicitly given). */

template<class A, class T, int Dim23, int Dim1, char i, char j, char k, int N0>
class Tensor3_dg_Expr<Tensor4_ddg_number_rhs_0<A,T,N0>,
  T,Dim23,Dim1,i,j,k>
{
  A &iter;
public:
  Tensor3_dg_Expr(A &a): iter(a) {}
  T & operator()(const int N1, const int N2, const int N3)
  {
    return iter(N0,N1,N2,N3);
  }
  T operator()(const int N1, const int N2, const int N3) const
  {
    return iter(N0,N1,N2,N3);
  }

  /* Various assignment operators.  I have to explicitly declare the
     second operator= because otherwise the compiler will generate its
     own and not use the template code. */

  template<class B, class U>
  const Tensor3_dg_Expr<Tensor4_ddg_number_rhs_0<A,T,N0>,
    T,Dim23,Dim1,i,j,k> &
  operator=(const Tensor3_dg_Expr<B,U,Dim23,Dim1,i,j,k> &result);

  const Tensor3_dg_Expr<Tensor4_ddg_number_rhs_0<A,T,N0>,
    T,Dim23,Dim1,i,j,k> &
  operator=(const Tensor3_dg_Expr<Tensor4_ddg_number_rhs_0<A,T,N0>,
	    T,Dim23,Dim1,i,j,k> &result);
};

#include "Tensor3_dg_Expr_equals.h"
