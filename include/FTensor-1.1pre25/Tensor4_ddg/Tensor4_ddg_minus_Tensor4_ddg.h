/* Subtracts Tensor4_ddg-Tensor4_ddg -> Tensor4_ddg */

/* A(i,j,k,l) - B(i,j,k,l) */

template<class A, class B, class T, class U, int Dim01, int Dim23,
  char i, char j, char k, char l>
class Tensor4_ddg_minus_Tensor4_ddg
{
  const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> iterA;
  const Tensor4_ddg_Expr<B,U,Dim01,Dim23,i,j,k,l> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2, const int N3,
			     const int N4) const
  {
    return iterA(N1,N2,N3,N4)-iterB(N1,N2,N3,N4);
  }

  Tensor4_ddg_minus_Tensor4_ddg(const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a, const Tensor4_ddg_Expr<B,U,Dim01,Dim23,i,j,k,l> &b):
    iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim01, int Dim23,
  char i, char j, char k, char l>
inline const Tensor4_ddg_Expr
<const Tensor4_ddg_minus_Tensor4_ddg<A,B,T,U,Dim01,Dim23,i,j,k,l>,
  typename promote<T,U>::V,Dim01,Dim23,i,j,k,l>
operator-(const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a,
	  const Tensor4_ddg_Expr<B,U,Dim01,Dim23,i,j,k,l> &b)
{
  typedef const Tensor4_ddg_minus_Tensor4_ddg<A,B,T,U,Dim01,Dim23,i,j,k,l>
    TensorExpr;
  return Tensor4_ddg_Expr<TensorExpr,typename promote<T,U>::V,Dim01,Dim23,i,j,k,l>
    (TensorExpr(a,b));
}
