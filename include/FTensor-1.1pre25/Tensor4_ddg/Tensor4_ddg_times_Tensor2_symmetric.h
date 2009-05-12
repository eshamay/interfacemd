/* This file has all of the declarations for expressions like
   Tensor4_ddg*Tensor2_symmetric and Tensor2_symmetric*Tensor4_ddg,
   yielding a Tensor2 or Tensor2_symmetric. */

/* A(i,j,k,l)*B(j,l) */

template<class A, class B, class T, class U, int Dim,
  char i, char j, char k, char l>
class Tensor4_ddg_times_Tensor2_symmetric_13
{
  const Tensor4_ddg_Expr<A,T,Dim,Dim,i,j,k,l> iterA;
  const Tensor2_symmetric_Expr<B,U,Dim,j,l> iterB;

  template<int Current_Dim0, int Current_Dim1>
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<Current_Dim0> &ND0,
		       const Number<Current_Dim1> &ND1)  const
  {
    return iterA(Current_Dim0-1,N1,Current_Dim1-1,N2)
      *iterB(Current_Dim0-1,Current_Dim1-1)
      + eval(N1,N2,Number<Current_Dim0-1>(),Number<Current_Dim1>());
  }
  template<int Current_Dim1>
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<1> &ND0,
		       const Number<Current_Dim1> &ND1)  const
  {
    return iterA(0,N1,Current_Dim1-1,N2)*iterB(0,Current_Dim1-1)
      + eval(N1,N2,Number<Dim>(),Number<Current_Dim1-1>());
  }
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<1> &ND0,
		       const Number<1> &ND1)  const
  {
    return iterA(0,N1,0,N2)*iterB(0,0);
  }
public:
  Tensor4_ddg_times_Tensor2_symmetric_13
  (const Tensor4_ddg_Expr<A,T,Dim,Dim,i,j,k,l> &a,
   const Tensor2_symmetric_Expr<B,U,Dim,j,l> &b): iterA(a), iterB(b) {}
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return eval(N1,N2,Number<Dim>(),Number<Dim>());
  }
};

template<class A, class B, class T, class U, int Dim,
  char i, char j, char k, char l>
inline const Tensor2_Expr
<const Tensor4_ddg_times_Tensor2_symmetric_13<A,B,T,U,Dim,i,j,k,l>,
  typename promote<T,U>::V,Dim,Dim,i,k>
operator*(const Tensor4_ddg_Expr<A,T,Dim,Dim,i,j,k,l> &a,
	  const Tensor2_symmetric_Expr<B,U,Dim,j,l> &b)
{
  typedef const Tensor4_ddg_times_Tensor2_symmetric_13<A,B,T,U,Dim,i,j,k,l>
    TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim,Dim,i,k>(TensorExpr(a,b));
}

/* B(j,l)*A(i,j,k,l) */

template<class A, class B, class T, class U, int Dim,
  char i, char j, char k, char l>
inline const Tensor2_Expr
<const Tensor4_ddg_times_Tensor2_symmetric_13<A,B,T,U,Dim,i,j,k,l>,
  typename promote<T,U>::V,Dim,Dim,i,k>
operator*(const Tensor2_symmetric_Expr<B,U,Dim,j,l> &b,
	  const Tensor4_ddg_Expr<A,T,Dim,Dim,i,j,k,l> &a)
{
  typedef const Tensor4_ddg_times_Tensor2_symmetric_13<A,B,T,U,Dim,i,j,k,l>
    TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim,Dim,i,k>(TensorExpr(a,b));
}

/* A(i,j,k,l)*B(i,j) */

template<class A, class B, class T, class U, int Dim01, int Dim23,
  char i, char j, char k, char l>
class Tensor4_ddg_times_Tensor2_symmetric_01
{
  const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> iterA;
  const Tensor2_symmetric_Expr<B,U,Dim01,i,j> iterB;

  template<int Current_Dim0, int Current_Dim1>
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<Current_Dim0> &ND0,
		       const Number<Current_Dim1> &ND1)  const
  {
    return iterA(Current_Dim0-1,Current_Dim1-1,N1,N2)
      *iterB(Current_Dim0-1,Current_Dim1-1)
      + eval(N1,N2,Number<Current_Dim0-1>(),Number<Current_Dim1>());
  }
  template<int Current_Dim1>
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<1> &ND0,
		       const Number<Current_Dim1> &ND1)  const
  {
    return iterA(0,Current_Dim1-1,N1,N2)*iterB(0,Current_Dim1-1)
      + eval(N1,N2,Number<Dim01>(),Number<Current_Dim1-1>());
  }
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<1> &ND0,
		       const Number<1> &ND1)  const
  {
    return iterA(0,0,N1,N2)*iterB(0,0);
  }
public:
  Tensor4_ddg_times_Tensor2_symmetric_01
  (const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a,
   const Tensor2_symmetric_Expr<B,U,Dim01,i,j> &b): iterA(a), iterB(b) {}
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return eval(N1,N2,Number<Dim01>(),Number<Dim01>());
  }
};

template<class A, class B, class T, class U, int Dim01, int Dim23,
  char i, char j, char k, char l>
inline const Tensor2_symmetric_Expr
<const Tensor4_ddg_times_Tensor2_symmetric_01<A,B,T,U,Dim01,Dim23,i,j,k,l>,
  typename promote<T,U>::V,Dim23,k,l>
operator*(const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a,
	  const Tensor2_symmetric_Expr<B,U,Dim01,i,j> &b)
{
  typedef const Tensor4_ddg_times_Tensor2_symmetric_01
    <A,B,T,U,Dim01,Dim23,i,j,k,l> TensorExpr;
  return Tensor2_symmetric_Expr<TensorExpr,typename promote<T,U>::V,Dim23,k,l>
    (TensorExpr(a,b));
}

/* B(i,j)*A(i,j,k,l) */

template<class A, class B, class T, class U, int Dim01, int Dim23,
  char i, char j, char k, char l>
inline const Tensor2_symmetric_Expr
<const Tensor4_ddg_times_Tensor2_symmetric_01<A,B,T,U,Dim01,Dim23,i,j,k,l>,
  typename promote<T,U>::V,Dim23,k,l>
operator*(const Tensor2_symmetric_Expr<B,U,Dim01,i,j> &b,
	  const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a)
{
  typedef const Tensor4_ddg_times_Tensor2_symmetric_01
    <A,B,T,U,Dim01,Dim23,i,j,k,l> TensorExpr;
  return Tensor2_symmetric_Expr<TensorExpr,typename promote<T,U>::V,Dim23,k,l>
    (TensorExpr(a,b));
}

/* A(i,j,k,l)*B(k,l) */

template<class A, class B, class T, class U, int Dim01, int Dim23,
  char i, char j, char k, char l>
class Tensor4_ddg_times_Tensor2_symmetric_23
{
  const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> iterA;
  const Tensor2_symmetric_Expr<B,U,Dim23,k,l> iterB;

  template<int Current_Dim0, int Current_Dim1>
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<Current_Dim0> &ND0,
		       const Number<Current_Dim1> &ND1)  const
  {
    return iterA(N1,N2,Current_Dim0-1,Current_Dim1-1)
      *iterB(Current_Dim0-1,Current_Dim1-1)
      + eval(N1,N2,Number<Current_Dim0-1>(),Number<Current_Dim1>());
  }
  template<int Current_Dim1>
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<1> &ND0,
		       const Number<Current_Dim1> &ND1)  const
  {
    return iterA(N1,N2,0,Current_Dim1-1)*iterB(0,Current_Dim1-1)
      + eval(N1,N2,Number<Dim23>(),Number<Current_Dim1-1>());
  }
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<1> &ND0,
		       const Number<1> &ND1)  const
  {
    return iterA(N1,N2,0,0)*iterB(0,0);
  }
public:
  Tensor4_ddg_times_Tensor2_symmetric_23
  (const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a,
   const Tensor2_symmetric_Expr<B,U,Dim23,k,l> &b): iterA(a), iterB(b) {}
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return eval(N1,N2,Number<Dim23>(),Number<Dim23>());
  }
};

template<class A, class B, class T, class U, int Dim01, int Dim23,
  char i, char j, char k, char l>
inline const Tensor2_symmetric_Expr
<const Tensor4_ddg_times_Tensor2_symmetric_23<A,B,T,U,Dim01,Dim23,i,j,k,l>,
  typename promote<T,U>::V,Dim01,i,j>
operator*(const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a,
	  const Tensor2_symmetric_Expr<B,U,Dim23,k,l> &b)
{
  typedef const Tensor4_ddg_times_Tensor2_symmetric_23
    <A,B,T,U,Dim01,Dim23,i,j,k,l> TensorExpr;
  return Tensor2_symmetric_Expr<TensorExpr,typename promote<T,U>::V,Dim01,i,j>
    (TensorExpr(a,b));
}

/* B(k,l)*A(i,j,k,l) */

template<class A, class B, class T, class U, int Dim01, int Dim23,
  char i, char j, char k, char l>
inline const Tensor2_symmetric_Expr
<const Tensor4_ddg_times_Tensor2_symmetric_23<A,B,T,U,Dim01,Dim23,i,j,k,l>,
  typename promote<T,U>::V,Dim01,i,j>
operator*(const Tensor2_symmetric_Expr<B,U,Dim23,k,l> &b,
	  const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a)
{
  typedef const Tensor4_ddg_times_Tensor2_symmetric_23
    <A,B,T,U,Dim01,Dim23,i,j,k,l> TensorExpr;
  return Tensor2_symmetric_Expr<TensorExpr,typename promote<T,U>::V,Dim01,i,j>
    (TensorExpr(a,b));
}

