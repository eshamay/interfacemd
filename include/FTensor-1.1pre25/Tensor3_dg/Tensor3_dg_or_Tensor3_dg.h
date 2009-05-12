/* Adds two Tensor3_dg's to make a Tensor3_dg, but with different
   symmetries. */

/* A(i,j,k)+B(i,k,j) -> Tensor3_dg(j,k,i) */

template<class A, class B, class T, class U, int Dim, char i, char j, char k>
class Tensor3_dg_or_Tensor3_dg_12
{
  const Tensor3_dg_Expr<A,T,Dim,Dim,i,j,k> iterA;
  const Tensor3_dg_Expr<B,U,Dim,Dim,i,k,j> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2, const int N3) const
  {
    return iterA(N3,N1,N2)+iterB(N3,N2,N1);
  }

  Tensor3_dg_or_Tensor3_dg_12(const Tensor3_dg_Expr<A,T,Dim,Dim,i,j,k> &a,
			      const Tensor3_dg_Expr<B,U,Dim,Dim,i,k,j> &b)
    : iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim, char i, char j, char k>
inline const Tensor3_dg_Expr<const Tensor3_dg_or_Tensor3_dg_12<A,B,T,U,Dim,i,j,k>,
  typename promote<T,U>::V,Dim,Dim,j,k,i>
operator||(const Tensor3_dg_Expr<A,T,Dim,Dim,i,j,k> &a,
	   const Tensor3_dg_Expr<B,U,Dim,Dim,i,k,j> &b)
{
  typedef const Tensor3_dg_or_Tensor3_dg_12<A,B,T,U,Dim,i,j,k> TensorExpr;
  return Tensor3_dg_Expr<TensorExpr,typename promote<T,U>::V,Dim,Dim,j,k,i>
    (TensorExpr(a,b));
}

/* A(j,i,k)+B(k,i,j) -> Tensor3_dg(j,k,i) */

template<class A, class B, class T, class U, int Dim, char i, char j, char k>
class Tensor3_dg_or_Tensor3_dg_02
{
  const Tensor3_dg_Expr<A,T,Dim,Dim,j,i,k> iterA;
  const Tensor3_dg_Expr<B,U,Dim,Dim,k,i,j> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2, const int N3) const
  {
    return iterA(N1,N3,N2)+iterB(N2,N3,N1);
  }

  Tensor3_dg_or_Tensor3_dg_02(const Tensor3_dg_Expr<A,T,Dim,Dim,j,i,k> &a,
			      const Tensor3_dg_Expr<B,U,Dim,Dim,k,i,j> &b)
    : iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim, char i, char j, char k>
inline const Tensor3_dg_Expr<const Tensor3_dg_or_Tensor3_dg_02<A,B,T,U,Dim,i,j,k>,
  typename promote<T,U>::V,Dim,Dim,j,k,i>
operator||(const Tensor3_dg_Expr<A,T,Dim,Dim,j,i,k> &a,
	   const Tensor3_dg_Expr<B,U,Dim,Dim,k,i,j> &b)
{
  typedef const Tensor3_dg_or_Tensor3_dg_02<A,B,T,U,Dim,i,j,k> TensorExpr;
  return Tensor3_dg_Expr<TensorExpr,typename promote<T,U>::V,Dim,Dim,j,k,i>
    (TensorExpr(a,b));
}



