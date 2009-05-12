/* Multiply a Tensor2_symmetric and a Tensor2_symmetric together but
   don't contract, yielding a Tensor2_symmetric. */

/* A(i,j) & B(i,j) -> Tensor2_symmetric */

template<class A, class B, class T, class U, int Dim, char i, char j>
class Tensor2_symmetric_and_Tensor2_symmetric_01
{
  const Tensor2_symmetric_Expr<A,T,Dim,i,j> iterA;
  const Tensor2_symmetric_Expr<B,U,Dim,i,j> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return iterA(N1,N2)*iterB(N1,N2);
  }

  Tensor2_symmetric_and_Tensor2_symmetric_01
  (const Tensor2_symmetric_Expr<A,T,Dim,i,j> &a,
   const Tensor2_symmetric_Expr<B,U,Dim,i,j> &b): iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim, char i, char j>
inline const Tensor2_symmetric_Expr
<const Tensor2_symmetric_and_Tensor2_symmetric_01<A,B,T,U,Dim,i,j>,
  typename promote<T,U>::V,Dim,i,j>
operator&(const Tensor2_symmetric_Expr<A,T,Dim,i,j> &a,
	  const Tensor2_symmetric_Expr<B,U,Dim,i,j> &b)
{
  typedef const Tensor2_symmetric_and_Tensor2_symmetric_01<A,B,T,U,Dim,i,j>
    TensorExpr;
  return Tensor2_symmetric_Expr<TensorExpr,typename promote<T,U>::V,Dim,i,j>
    (TensorExpr(a,b));
}

/* A(i,j) & B(j,i) -> Tensor2_symmetric */

template<class A, class B, class T, class U, int Dim, char i, char j>
class Tensor2_symmetric_and_Tensor2_symmetric_10
{
  const Tensor2_symmetric_Expr<A,T,Dim,i,j> iterA;
  const Tensor2_symmetric_Expr<B,U,Dim,j,i> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return iterA(N1,N2)*iterB(N2,N1);
  }

  Tensor2_symmetric_and_Tensor2_symmetric_10
  (const Tensor2_symmetric_Expr<A,T,Dim,i,j> &a,
   const Tensor2_symmetric_Expr<B,U,Dim,j,i> &b): iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim, char i, char j>
inline const Tensor2_symmetric_Expr<const Tensor2_symmetric_and_Tensor2_symmetric_10
<A,B,T,U,Dim,i,j>,typename promote<T,U>::V,Dim,i,j>
operator&(const Tensor2_symmetric_Expr<A,T,Dim,i,j> &a,
	  const Tensor2_symmetric_Expr<B,U,Dim,j,i> &b)
{
  typedef const Tensor2_symmetric_and_Tensor2_symmetric_10<A,B,T,U,Dim,i,j>
    TensorExpr;
  return Tensor2_symmetric_Expr<TensorExpr,typename promote<T,U>::V,Dim,i,j>
    (TensorExpr(a,b));
}

