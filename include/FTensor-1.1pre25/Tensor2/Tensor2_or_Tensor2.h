/* Creates a Tensor2_symmetric expression by adding two Tensor2's
   together. */

/* This version is for expressions like A(i,j)||B(i,j) */

template<class A, class B, class T, class U, int Dim, char i, char j>
class Tensor2_or_Tensor2
{
  const Tensor2_Expr<A,T,Dim,Dim,i,j> iterA;
  const Tensor2_Expr<B,U,Dim,Dim,i,j> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return iterA(N1,N2)+iterB(N1,N2);
  }

  Tensor2_or_Tensor2(const Tensor2_Expr<A,T,Dim,Dim,i,j> &a,
		     const Tensor2_Expr<B,U,Dim,Dim,i,j> &b):
    iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim, char i, char j>
inline const Tensor2_symmetric_Expr<const Tensor2_or_Tensor2<A,B,T,U,Dim,i,j>,
  typename promote<T,U>::V,Dim,i,j>
operator||(const Tensor2_Expr<A,T,Dim,Dim,i,j> &a,
	   const Tensor2_Expr<B,U,Dim,Dim,i,j> &b)
{
  typedef const Tensor2_or_Tensor2<A,B,T,U,Dim,i,j> TensorExpr;
  return Tensor2_symmetric_Expr<TensorExpr,typename promote<T,U>::V,Dim,i,j>
    (TensorExpr(a,b));
}

/* This version is for when the indices are switched as in
   A(i,j)||B(j,i) (probably the most common case). */

template<class A, class B, class T, class U, int Dim, char i, char j>
class Tensor2_or_Tensor2_switched
{
  const Tensor2_Expr<A,T,Dim,Dim,i,j> iterA;
  const Tensor2_Expr<B,U,Dim,Dim,j,i> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return iterA(N1,N2)+iterB(N2,N1);
  }

  Tensor2_or_Tensor2_switched(const Tensor2_Expr<A,T,Dim,Dim,i,j> &a,
			      const Tensor2_Expr<B,U,Dim,Dim,j,i> &b):
    iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim, char i, char j>
inline const Tensor2_symmetric_Expr
<const Tensor2_or_Tensor2_switched<A,B,T,U,Dim,i,j>,typename promote<T,U>::V,Dim,i,j>
operator||(const Tensor2_Expr<A,T,Dim,Dim,i,j> &a,
	   const Tensor2_Expr<B,U,Dim,Dim,j,i> &b)
{
  typedef const Tensor2_or_Tensor2_switched<A,B,T,U,Dim,i,j> TensorExpr;
  return Tensor2_symmetric_Expr<TensorExpr,typename promote<T,U>::V,Dim,i,j>
    (TensorExpr(a,b));
}
