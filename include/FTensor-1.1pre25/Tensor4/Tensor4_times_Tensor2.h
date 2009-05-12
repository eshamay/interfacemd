/* This file has all of the declarations for expressions like
   Tensor4*Tensor2 and Tensor2*Tensor4, yielding a
   Tensor2. */

/* A(i,j,k,l)*B(k,l) */

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  int Dim3, char i, char j, char k, char l>
class Tensor4_times_Tensor2_23
{
  const Tensor4_Expr<A,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> iterA;
  const Tensor2_Expr<B,U,Dim2,Dim3,k,l> iterB;

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
      + eval(N1,N2,Number<Dim2>(),Number<Current_Dim1-1>());
  }
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<1> &ND0,
		       const Number<1> &ND1)  const
  {
    return iterA(N1,N2,0,0)*iterB(0,0);
  }
public:
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return eval(N1,N2,Number<Dim2>(),Number<Dim3>());
  }

  Tensor4_times_Tensor2_23
  (const Tensor4_Expr<A,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &a,
   const Tensor2_Expr<B,U,Dim2,Dim3,k,l> &b): iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  int Dim3, char i, char j, char k, char l>
inline Tensor2_Expr
<const Tensor4_times_Tensor2_23<A,B,T,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l>,
  typename promote<T,U>::V,Dim0,Dim1,i,j>
operator*(const Tensor4_Expr<A,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &a,
	  const Tensor2_Expr<B,U,Dim2,Dim3,k,l> &b)
{
  typedef const Tensor4_times_Tensor2_23
    <A,B,T,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim0,Dim1,i,j>
    (TensorExpr(a,b));
}

/* B(k,l)*A(i,j,k,l) */

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  int Dim3, char i, char j, char k, char l>
inline Tensor2_Expr
<const Tensor4_times_Tensor2_23<A,B,T,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l>,
  typename promote<T,U>::V,Dim0,Dim1,i,j>
operator*(const Tensor2_Expr<B,U,Dim2,Dim3,k,l> &b,
	  const Tensor4_Expr<A,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &a)
{
  typedef const Tensor4_times_Tensor2_23
    <A,B,T,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim0,Dim1,i,j>
    (TensorExpr(a,b));
}

/* A(i,j,k,l)*B(l,k) */

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  int Dim3, char i, char j, char k, char l>
class Tensor4_times_Tensor2_32
{
  const Tensor4_Expr<A,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> iterA;
  const Tensor2_Expr<B,U,Dim3,Dim2,l,k> iterB;

  template<int Current_Dim0, int Current_Dim1>
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<Current_Dim0> &ND0,
		       const Number<Current_Dim1> &ND1)  const
  {
    return iterA(N1,N2,Current_Dim0-1,Current_Dim1-1)
      *iterB(Current_Dim1-1,Current_Dim0-1)
      + eval(N1,N2,Number<Current_Dim0-1>(),Number<Current_Dim1>());
  }
  template<int Current_Dim1>
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<1> &ND0,
		       const Number<Current_Dim1> &ND1)  const
  {
    return iterA(N1,N2,0,Current_Dim1-1)*iterB(Current_Dim1-1,0)
      + eval(N1,N2,Number<Dim2>(),Number<Current_Dim1-1>());
  }
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<1> &ND0,
		       const Number<1> &ND1)  const
  {
    return iterA(N1,N2,0,0)*iterB(0,0);
  }
public:
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return eval(N1,N2,Number<Dim2>(),Number<Dim3>());
  }

  Tensor4_times_Tensor2_32
  (const Tensor4_Expr<A,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &a,
   const Tensor2_Expr<B,U,Dim3,Dim2,l,k> &b): iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  int Dim3, char i, char j, char k, char l>
inline Tensor2_Expr
<const Tensor4_times_Tensor2_32<A,B,T,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l>,
  typename promote<T,U>::V,Dim0,Dim1,i,j>
operator*(const Tensor4_Expr<A,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &a,
	  const Tensor2_Expr<B,U,Dim3,Dim2,l,k> &b)
{
  typedef const Tensor4_times_Tensor2_32
    <A,B,T,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim0,Dim1,i,j>
    (TensorExpr(a,b));
}

/* B(l,k)*A(i,j,k,l) */

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  int Dim3, char i, char j, char k, char l>
inline Tensor2_Expr
<const Tensor4_times_Tensor2_32<A,B,T,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l>,
  typename promote<T,U>::V,Dim0,Dim1,i,j>
operator*(const Tensor2_Expr<B,U,Dim3,Dim2,l,k> &b,
	  const Tensor4_Expr<A,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &a)
{
  typedef const Tensor4_times_Tensor2_32
    <A,B,T,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim0,Dim1,i,j>
    (TensorExpr(a,b));
}

/* A(i,j,k,l)*B(i,l) */

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  int Dim3, char i, char j, char k, char l>
class Tensor4_times_Tensor2_03
{
  const Tensor4_Expr<A,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> iterA;
  const Tensor2_Expr<B,U,Dim0,Dim3,i,l> iterB;

  template<int Current_Dim0, int Current_Dim1>
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<Current_Dim0> &ND0,
		       const Number<Current_Dim1> &ND1)  const
  {
    return iterA(Current_Dim0-1,N1,N2,Current_Dim1-1)
      *iterB(Current_Dim0-1,Current_Dim1-1)
      + eval(N1,N2,Number<Current_Dim0-1>(),Number<Current_Dim1>());
  }
  template<int Current_Dim1>
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<1> &ND0,
		       const Number<Current_Dim1> &ND1)  const
  {
    return iterA(0,N1,N2,Current_Dim1-1)*iterB(0,Current_Dim1-1)
      + eval(N1,N2,Number<Dim0>(),Number<Current_Dim1-1>());
  }
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<1> &ND0,
		       const Number<1> &ND1)  const
  {
    return iterA(0,N1,N2,0)*iterB(0,0);
  }
public:
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return eval(N1,N2,Number<Dim0>(),Number<Dim3>());
  }

  Tensor4_times_Tensor2_03
  (const Tensor4_Expr<A,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &a,
   const Tensor2_Expr<B,U,Dim0,Dim3,i,l> &b): iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  int Dim3, char i, char j, char k, char l>
inline Tensor2_Expr
<const Tensor4_times_Tensor2_03<A,B,T,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l>,
  typename promote<T,U>::V,Dim1,Dim2,j,k>
operator*(const Tensor4_Expr<A,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &a,
	  const Tensor2_Expr<B,U,Dim0,Dim3,i,l> &b)
{
  typedef const Tensor4_times_Tensor2_03
    <A,B,T,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim1,Dim2,j,k>
    (TensorExpr(a,b));
}

/* B(i,l)*A(i,j,k,l) */

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  int Dim3, char i, char j, char k, char l>
inline Tensor2_Expr
<const Tensor4_times_Tensor2_03<A,B,T,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l>,
  typename promote<T,U>::V,Dim1,Dim2,j,k>
operator*(const Tensor2_Expr<B,U,Dim0,Dim3,i,l> &b,
	  const Tensor4_Expr<A,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &a)
{
  typedef const Tensor4_times_Tensor2_03
    <A,B,T,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim1,Dim2,j,k>
    (TensorExpr(a,b));
}

/* A(i,j,k,l)*B(l,i) */

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  int Dim3, char i, char j, char k, char l>
class Tensor4_times_Tensor2_30
{
  const Tensor4_Expr<A,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> iterA;
  const Tensor2_Expr<B,U,Dim3,Dim0,l,i> iterB;

  template<int Current_Dim0, int Current_Dim1>
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<Current_Dim0> &ND0,
		       const Number<Current_Dim1> &ND1)  const
  {
    return iterA(Current_Dim0-1,N1,N2,Current_Dim1-1)
      *iterB(Current_Dim1-1,Current_Dim0-1)
      + eval(N1,N2,Number<Current_Dim0-1>(),Number<Current_Dim1>());
  }
  template<int Current_Dim1>
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<1> &ND0,
		       const Number<Current_Dim1> &ND1)  const
  {
    return iterA(0,N1,N2,Current_Dim1-1)*iterB(Current_Dim1-1,0)
      + eval(N1,N2,Number<Dim0>(),Number<Current_Dim1-1>());
  }
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<1> &ND0,
		       const Number<1> &ND1)  const
  {
    return iterA(0,N1,N2,0)*iterB(0,0);
  }
public:
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return eval(N1,N2,Number<Dim0>(),Number<Dim3>());
  }

  Tensor4_times_Tensor2_30
  (const Tensor4_Expr<A,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &a,
   const Tensor2_Expr<B,U,Dim3,Dim0,l,i> &b): iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  int Dim3, char i, char j, char k, char l>
inline Tensor2_Expr
<const Tensor4_times_Tensor2_30<A,B,T,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l>,
  typename promote<T,U>::V,Dim1,Dim2,j,k>
operator*(const Tensor4_Expr<A,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &a,
	  const Tensor2_Expr<B,U,Dim3,Dim0,l,i> &b)
{
  typedef const Tensor4_times_Tensor2_30
    <A,B,T,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim1,Dim2,j,k>
    (TensorExpr(a,b));
}

/* B(l,i)*A(i,j,k,l) */

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  int Dim3, char i, char j, char k, char l>
inline Tensor2_Expr
<const Tensor4_times_Tensor2_30<A,B,T,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l>,
  typename promote<T,U>::V,Dim1,Dim2,j,k>
operator*(const Tensor2_Expr<B,U,Dim3,Dim0,l,i> &b,
	  const Tensor4_Expr<A,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &a)
{
  typedef const Tensor4_times_Tensor2_30
    <A,B,T,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim1,Dim2,j,k>
    (TensorExpr(a,b));
}

