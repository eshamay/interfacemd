/* Adds a Tensor4_Riemann to a Tensor4_Riemann, yielding a
   Tensor4_Riemann. */

template<class A, class B, class T, class U, int Dim,
  char i, char j, char k, char l>
class Tensor4_Riemann_plus_Tensor4_Riemann
{
  const Tensor4_Riemann_Expr<A,T,Dim,i,j,k,l> iterA;
  const Tensor4_Riemann_Expr<B,U,Dim,i,j,k,l> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2, const int N3,
			     const int N4) const
  {
    return iterA(N1,N2,N3,N4)+iterB(N1,N2,N3,N4);
  }

  Tensor4_Riemann_plus_Tensor4_Riemann
  (const Tensor4_Riemann_Expr<A,T,Dim,i,j,k,l> &a,
   const Tensor4_Riemann_Expr<B,U,Dim,i,j,k,l> &b): iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim,
  char i, char j, char k, char l>
inline const Tensor4_Riemann_Expr
<const Tensor4_Riemann_plus_Tensor4_Riemann<A,B,T,U,Dim,i,j,k,l>,
  typename promote<T,U>::V,Dim,i,j,k,l>
operator+(const Tensor4_Riemann_Expr<A,T,Dim,i,j,k,l> &a,
	  const Tensor4_Riemann_Expr<B,U,Dim,i,j,k,l> &b)
{
  typedef const Tensor4_Riemann_plus_Tensor4_Riemann<A,B,T,U,Dim,i,j,k,l>
    TensorExpr;
  return Tensor4_Riemann_Expr<TensorExpr,typename promote<T,U>::V,Dim,i,j,k,l>
    (TensorExpr(a,b));
}
