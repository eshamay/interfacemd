/* Declares expressions of Tensor4_ddg || Tensor4_ddg.  This adds them in
   a different way, but still ending up with a Tensor4_ddg. */

/* A(i,j,k,l)+B(i,l,k,j) -> Tensor4_ddg */

template<class A, class B, class T, class U, int Dim,
  char i, char j, char k, char l>
class Tensor4_ddg_or_Tensor4_ddg
{
  const Tensor4_ddg_Expr<A,T,Dim,Dim,i,j,k,l> iterA;
  const Tensor4_ddg_Expr<B,U,Dim,Dim,i,l,k,j> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2, const int N3,
			     const int N4) const
  {
    return iterA(N1,N3,N2,N4)+iterB(N1,N4,N2,N3);
  }
  Tensor4_ddg_or_Tensor4_ddg(const Tensor4_ddg_Expr<A,T,Dim,Dim,i,j,k,l> &a,
			     const Tensor4_ddg_Expr<B,U,Dim,Dim,i,l,k,j> &b)
    : iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim,
  char i, char j, char k, char l>
inline const Tensor4_ddg_Expr
<const Tensor4_ddg_or_Tensor4_ddg<A,B,T,U,Dim,i,j,k,l>,
  typename promote<T,U>::V,Dim,Dim,i,k,j,l>
operator||(const Tensor4_ddg_Expr<A,T,Dim,Dim,i,j,k,l> &a,
	   const Tensor4_ddg_Expr<B,U,Dim,Dim,i,l,k,j> &b)
{
  typedef const Tensor4_ddg_or_Tensor4_ddg<A,B,T,U,Dim,i,j,k,l> TensorExpr;
  return Tensor4_ddg_Expr<TensorExpr,typename promote<T,U>::V,Dim,Dim,i,k,j,l>
    (TensorExpr(a,b));
}
