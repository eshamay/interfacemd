/* Adds a Tensor2_symmetric to a Tensor2, yielding a Tensor2. */

/* A(i,j)+B(i,j), A is symmetric, B is not. */

template<class A, class B, class T, class U, int Dim, char i, char j>
class Tensor2_symmetric_plus_Tensor2_01
{
  const Tensor2_symmetric_Expr<A,T,Dim,i,j> iterA;
  const Tensor2_Expr<B,U,Dim,Dim,i,j> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return iterA(N1,N2)+iterB(N1,N2);
  }

  Tensor2_symmetric_plus_Tensor2_01
  (const Tensor2_symmetric_Expr<A,T,Dim,i,j> &a,
   const Tensor2_Expr<B,U,Dim,Dim,i,j> &b): iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim, char i, char j>
inline const Tensor2_Expr<const Tensor2_symmetric_plus_Tensor2_01<A,B,T,U,Dim,i,j>,
  typename promote<T,U>::V,Dim,Dim,i,j>
operator+(const Tensor2_symmetric_Expr<A,T,Dim,i,j> &a,
	  const Tensor2_Expr<B,U,Dim,Dim,i,j> &b)
{
  typedef const Tensor2_symmetric_plus_Tensor2_01<A,B,T,U,Dim,i,j> TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim,Dim,i,j>(TensorExpr(a,b));
}

/* B(i,j)+A(i,j), A is symmetric, B is not. */

template<class A, class B, class T, class U, int Dim, char i, char j>
inline const Tensor2_Expr<const Tensor2_symmetric_plus_Tensor2_01<A,B,T,U,Dim,i,j>,
  typename promote<T,U>::V,Dim,Dim,i,j>
operator+(const Tensor2_Expr<B,U,Dim,Dim,i,j> &b,
	  const Tensor2_symmetric_Expr<A,T,Dim,i,j> &a)
{
  typedef const Tensor2_symmetric_plus_Tensor2_01<A,B,T,U,Dim,i,j> TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim,Dim,i,j>(TensorExpr(a,b));
}

/* A(i,j)+B(j,i), A is symmetric, B is not. */

template<class A, class B, class T, class U, int Dim, char i, char j>
class Tensor2_symmetric_plus_Tensor2_10
{
  const Tensor2_symmetric_Expr<A,T,Dim,i,j> iterA;
  const Tensor2_Expr<B,U,Dim,Dim,j,i> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return iterA(N1,N2)+iterB(N2,N1);
  }

  Tensor2_symmetric_plus_Tensor2_10
  (const Tensor2_symmetric_Expr<A,T,Dim,i,j> &a,
   const Tensor2_Expr<B,U,Dim,Dim,j,i> &b): iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim, char i, char j>
inline const Tensor2_Expr<const Tensor2_symmetric_plus_Tensor2_10<A,B,T,U,Dim,i,j>,
  typename promote<T,U>::V,Dim,Dim,i,j>
operator+(const Tensor2_symmetric_Expr<A,T,Dim,i,j> &a,
	  const Tensor2_Expr<B,U,Dim,Dim,j,i> &b)
{
  typedef const Tensor2_symmetric_plus_Tensor2_10<A,B,T,U,Dim,i,j> TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim,Dim,i,j>(TensorExpr(a,b));
}

/* B(j,i)+A(i,j), A is symmetric, B is not. */

template<class A, class B, class T, class U, int Dim, char i, char j>
inline const Tensor2_Expr<const Tensor2_symmetric_plus_Tensor2_10<A,B,T,U,Dim,i,j>,
  typename promote<T,U>::V,Dim,Dim,i,j>
operator+(const Tensor2_Expr<B,U,Dim,Dim,j,i> &b,
	  const Tensor2_symmetric_Expr<A,T,Dim,i,j> &a)
{
  typedef const Tensor2_symmetric_plus_Tensor2_10<A,B,T,U,Dim,i,j> TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim,Dim,i,j>(TensorExpr(a,b));
}

