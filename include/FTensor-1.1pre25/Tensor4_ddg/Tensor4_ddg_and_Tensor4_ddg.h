/* Subtracts Tensor4_ddg-Tensor4_ddg -> Tensor4_Riemann */

/* A(i,j,k,l) - B(i,l,k,j) */

template<class A, class B, class T, class U, int Dim,
  char i, char j, char k, char l>
class Tensor4_ddg_and_Tensor4_ddg0321
{
  const Tensor4_ddg_Expr<A,T,Dim,Dim,i,j,k,l> iterA;
  const Tensor4_ddg_Expr<B,U,Dim,Dim,i,l,k,j> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2, const int N3,
			     const int N4) const
  {
    return iterA(N1,N2,N3,N4)-iterB(N1,N4,N3,N2);
  }

  Tensor4_ddg_and_Tensor4_ddg0321
  (const Tensor4_ddg_Expr<A,T,Dim,Dim,i,j,k,l> &a,
   const Tensor4_ddg_Expr<B,U,Dim,Dim,i,l,k,j> &b): iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim,
  char i, char j, char k, char l>
inline const Tensor4_Riemann_Expr
<const Tensor4_ddg_and_Tensor4_ddg0321<A,B,T,U,Dim,i,j,k,l>,
  typename promote<T,U>::V,Dim,i,j,k,l>
operator&&(const Tensor4_ddg_Expr<A,T,Dim,Dim,i,j,k,l> &a,
	  const Tensor4_ddg_Expr<B,U,Dim,Dim,i,l,k,j> &b)
{
  typedef const Tensor4_ddg_and_Tensor4_ddg0321<A,B,T,U,Dim,i,j,k,l> TensorExpr;
  return Tensor4_Riemann_Expr<TensorExpr,typename promote<T,U>::V,Dim,i,j,k,l>
    (TensorExpr(a,b));
}

/* A(i,k,l,j) - B(i,l,k,j) */

template<class A, class B, class T, class U, int Dim,
  char i, char j, char k, char l>
class Tensor4_ddg_and_Tensor4_ddg0213
{
  const Tensor4_ddg_Expr<A,T,Dim,Dim,i,k,l,j> iterA;
  const Tensor4_ddg_Expr<B,U,Dim,Dim,i,l,k,j> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2, const int N3,
			     const int N4) const
  {
    return iterA(N1,N3,N4,N2)-iterB(N1,N4,N3,N2);
  }

  Tensor4_ddg_and_Tensor4_ddg0213
  (const Tensor4_ddg_Expr<A,T,Dim,Dim,i,k,l,j> &a,
   const Tensor4_ddg_Expr<B,U,Dim,Dim,i,l,k,j> &b): iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim,
  char i, char j, char k, char l>
inline const Tensor4_Riemann_Expr
<const Tensor4_ddg_and_Tensor4_ddg0213<A,B,T,U,Dim,i,j,k,l>,
  typename promote<T,U>::V,Dim,i,j,k,l>
operator&&(const Tensor4_ddg_Expr<A,T,Dim,Dim,i,k,l,j> &a,
	  const Tensor4_ddg_Expr<B,U,Dim,Dim,i,l,k,j> &b)
{
  typedef const Tensor4_ddg_and_Tensor4_ddg0213<A,B,T,U,Dim,i,j,k,l>
    TensorExpr;
  return Tensor4_Riemann_Expr<TensorExpr,typename promote<T,U>::V,Dim,i,j,k,l>
    (TensorExpr(a,b));
}
