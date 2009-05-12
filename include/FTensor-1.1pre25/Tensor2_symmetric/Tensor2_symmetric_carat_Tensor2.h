/* Creates a Tensor2_symmetric expression by contracting a
   Tensor2_symmetric and a Tensor2 together. There are different
   versions, depending on where the contracting indices are located
   (i.e. whether it is A(i,j)^B(j,k) or A(i,j)^B(k,j)).  The classes
   are numbered to differentiate between these.  Thus, A(i,j)^B(j,k)
   has 10 appended to the name because I count from 0. */

/* A(i,j)*B(j,k) */

template<class A, class B, class T, class U, int Dim,
  char i, char j, char k>
class Tensor2_symmetric_carat_Tensor2_10
{
  const Tensor2_symmetric_Expr<A,T,Dim,i,j> iterA;
  const Tensor2_Expr<B,U,Dim,Dim,j,k> iterB;

  template<int Current_Dim>
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<Current_Dim> &N) const
  {
    return iterA(N1,Current_Dim-1)*iterB(Current_Dim-1,N2)
      + eval(N1,N2,Number<Current_Dim-1>());
  }
  typename promote<T,U>::V eval(const int N1, const int N2, const Number<1> &N) const
  {
    return iterA(N1,0)*iterB(0,N2);
  }
public:
  Tensor2_symmetric_carat_Tensor2_10
  (const Tensor2_symmetric_Expr<A,T,Dim,i,j> &a,
   const Tensor2_Expr<B,U,Dim,Dim,j,k> &b): iterA(a), iterB(b) {}
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return eval(N1,N2,Number<Dim>());
  }
};

template<class A, class B, class T, class U, int Dim,
  char i, char j, char k>
inline const Tensor2_symmetric_Expr
<const Tensor2_symmetric_carat_Tensor2_10<A,B,T,U,Dim,i,j,k>,
  typename promote<T,U>::V,Dim,i,k>
operator^(const Tensor2_symmetric_Expr<A,T,Dim,i,j> &a,
	  const Tensor2_Expr<B,U,Dim,Dim,j,k> &b)
{
  typedef const Tensor2_symmetric_carat_Tensor2_10<A,B,T,U,Dim,i,j,k>
    TensorExpr;
  return Tensor2_symmetric_Expr<TensorExpr,typename promote<T,U>::V,Dim,i,k>
    (TensorExpr(a,b));
}

/* B(j,k)*A(i,j) */

template<class A, class B, class T, class U, int Dim,
  char i, char j, char k>
inline const Tensor2_symmetric_Expr
<const Tensor2_symmetric_carat_Tensor2_10<A,B,T,U,Dim,i,j,k>,
  typename promote<T,U>::V,Dim,i,k>
operator^(const Tensor2_Expr<B,U,Dim,Dim,j,k> &b,
	  const Tensor2_symmetric_Expr<A,T,Dim,i,j> &a)
{
  typedef const Tensor2_symmetric_carat_Tensor2_10<A,B,T,U,Dim,i,j,k>
    TensorExpr;
  return Tensor2_symmetric_Expr<TensorExpr,typename promote<T,U>::V,Dim,i,k>
    (TensorExpr(a,b));
}







/* A(i,j)*B(k,j) */

template<class A, class B, class T, class U, int Dim,
  char i, char j, char k>
class Tensor2_symmetric_carat_Tensor2_11
{
  const Tensor2_symmetric_Expr<A,T,Dim,i,j> iterA;
  const Tensor2_Expr<B,U,Dim,Dim,k,j> iterB;

  template<int Current_Dim>
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<Current_Dim> &N) const
  {
    return iterA(N1,Current_Dim-1)*iterB(N2,Current_Dim-1)
      + eval(N1,N2,Number<Current_Dim-1>());
  }
  typename promote<T,U>::V eval(const int N1, const int N2, const Number<1> &N) const
  {
    return iterA(N1,0)*iterB(N2,0);
  }
public:
  Tensor2_symmetric_carat_Tensor2_11
  (const Tensor2_symmetric_Expr<A,T,Dim,i,j> &a,
   const Tensor2_Expr<B,U,Dim,Dim,k,j> &b): iterA(a), iterB(b) {}
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return eval(N1,N2,Number<Dim>());
  }
};

template<class A, class B, class T, class U, int Dim,
  char i, char j, char k>
inline const Tensor2_symmetric_Expr
<const Tensor2_symmetric_carat_Tensor2_11<A,B,T,U,Dim,i,j,k>,
  typename promote<T,U>::V,Dim,i,k>
operator^(const Tensor2_symmetric_Expr<A,T,Dim,i,j> &a,
	  const Tensor2_Expr<B,U,Dim,Dim,k,j> &b)
{
  typedef const Tensor2_symmetric_carat_Tensor2_11<A,B,T,U,Dim,i,j,k>
    TensorExpr;
  return Tensor2_symmetric_Expr<TensorExpr,typename promote<T,U>::V,Dim,i,k>
    (TensorExpr(a,b));
}

/* B(k,j)*A(i,j) */

template<class A, class B, class T, class U, int Dim,
  char i, char j, char k>
inline const Tensor2_symmetric_Expr
<const Tensor2_symmetric_carat_Tensor2_11<A,B,T,U,Dim,i,j,k>,
  typename promote<T,U>::V,Dim,i,k>
operator^(const Tensor2_Expr<B,U,Dim,Dim,k,j> &b,
	  const Tensor2_symmetric_Expr<A,T,Dim,i,j> &a)
{
  typedef const Tensor2_symmetric_carat_Tensor2_11<A,B,T,U,Dim,i,j,k>
    TensorExpr;
  return Tensor2_symmetric_Expr<TensorExpr,typename promote<T,U>::V,Dim,i,k>
    (TensorExpr(a,b));
}







/* A(j,i)*B(j,k) */

template<class A, class B, class T, class U, int Dim,
  char i, char j, char k>
class Tensor2_symmetric_carat_Tensor2_00
{
  const Tensor2_symmetric_Expr<A,T,Dim,j,i> iterA;
  const Tensor2_Expr<B,U,Dim,Dim,j,k> iterB;

  template<int Current_Dim>
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<Current_Dim> &N) const
  {
    return iterA(Current_Dim-1,N1)*iterB(Current_Dim-1,N2)
      + eval(N1,N2,Number<Current_Dim-1>());
  }
  typename promote<T,U>::V eval(const int N1, const int N2, const Number<1> &N) const
  {
    return iterA(0,N1)*iterB(0,N2);
  }
public:
  Tensor2_symmetric_carat_Tensor2_00
  (const Tensor2_symmetric_Expr<A,T,Dim,j,i> &a,
   const Tensor2_Expr<B,U,Dim,Dim,j,k> &b): iterA(a), iterB(b) {}
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return eval(N1,N2,Number<Dim>());
  }
};

template<class A, class B, class T, class U, int Dim,
  char i, char j, char k>
inline const Tensor2_symmetric_Expr
<const Tensor2_symmetric_carat_Tensor2_00<A,B,T,U,Dim,i,j,k>,
  typename promote<T,U>::V,Dim,i,k>
operator^(const Tensor2_symmetric_Expr<A,T,Dim,j,i> &a,
	  const Tensor2_Expr<B,U,Dim,Dim,j,k> &b)
{
  typedef const Tensor2_symmetric_carat_Tensor2_00<A,B,T,U,Dim,i,j,k>
    TensorExpr;
  return Tensor2_symmetric_Expr<TensorExpr,typename promote<T,U>::V,Dim,i,k>
    (TensorExpr(a,b));
}

/* B(j,k)*A(j,i) */

template<class A, class B, class T, class U, int Dim,
  char i, char j, char k>
inline const Tensor2_symmetric_Expr
<const Tensor2_symmetric_carat_Tensor2_00<A,B,T,U,Dim,i,j,k>,
  typename promote<T,U>::V,Dim,i,k>
operator^(const Tensor2_Expr<B,U,Dim,Dim,j,k> &b,
	  const Tensor2_symmetric_Expr<A,T,Dim,j,i> &a)
{
  typedef const Tensor2_symmetric_carat_Tensor2_00<A,B,T,U,Dim,i,j,k>
    TensorExpr;
  return Tensor2_symmetric_Expr<TensorExpr,typename promote<T,U>::V,Dim,i,k>
    (TensorExpr(a,b));
}






/* A(j,i)*B(k,j) */

template<class A, class B, class T, class U, int Dim,
  char i, char j, char k>
class Tensor2_symmetric_carat_Tensor2_01
{
  const Tensor2_symmetric_Expr<A,T,Dim,j,i> iterA;
  const Tensor2_Expr<B,U,Dim,Dim,k,j> iterB;

  template<int Current_Dim>
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<Current_Dim> &N) const
  {
    return iterA(Current_Dim-1,N1)*iterB(N2,Current_Dim-1)
      + eval(N1,N2,Number<Current_Dim-1>());
  }
  typename promote<T,U>::V eval(const int N1, const int N2, const Number<1> &N) const
  {
    return iterA(0,N1)*iterB(N2,0);
  }
public:
  Tensor2_symmetric_carat_Tensor2_01
  (const Tensor2_symmetric_Expr<A,T,Dim,j,i> &a,
   const Tensor2_Expr<B,U,Dim,Dim,k,j> &b): iterA(a), iterB(b) {}
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return eval(N1,N2,Number<Dim>());
  }
};

template<class A, class B, class T, class U, int Dim,
  char i, char j, char k>
inline const Tensor2_symmetric_Expr
<const Tensor2_symmetric_carat_Tensor2_01<A,B,T,U,Dim,i,j,k>,
  typename promote<T,U>::V,Dim,i,k>
operator^(const Tensor2_symmetric_Expr<A,T,Dim,j,i> &a,
	  const Tensor2_Expr<B,U,Dim,Dim,k,j> &b)
{
  typedef const Tensor2_symmetric_carat_Tensor2_01<A,B,T,U,Dim,i,j,k>
    TensorExpr;
  return Tensor2_symmetric_Expr<TensorExpr,typename promote<T,U>::V,Dim,i,k>
    (TensorExpr(a,b));
}

/* B(k,j)*A(j,i) */

template<class A, class B, class T, class U, int Dim,
  char i, char j, char k>
inline const Tensor2_symmetric_Expr
<const Tensor2_symmetric_carat_Tensor2_01<A,B,T,U,Dim,i,j,k>,
  typename promote<T,U>::V,Dim,i,k>
operator^(const Tensor2_Expr<B,U,Dim,Dim,k,j> &b,
	  const Tensor2_symmetric_Expr<A,T,Dim,j,i> &a)
{
  typedef const Tensor2_symmetric_carat_Tensor2_01<A,B,T,U,Dim,i,j,k>
    TensorExpr;
  return Tensor2_symmetric_Expr<TensorExpr,typename promote<T,U>::V,Dim,i,k>
    (TensorExpr(a,b));
}
