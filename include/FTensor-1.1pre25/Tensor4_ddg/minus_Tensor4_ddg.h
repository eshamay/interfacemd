/* Unary minus operator. */

template<class A, class T, int Dim01, int Dim23, char i, char j, char k, char l>
class minus_Tensor4_ddg
{
public:
  const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> iterA;
public:
  T operator()(const int N1, const int N2, const int N3, const int N4) const
  {
    return -iterA(N1,N2,N3,N4);
  }
  minus_Tensor4_ddg(const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a):
    iterA(a) {}
};

template<class A, class T, int Dim01, int Dim23, char i, char j, char k, char l>
inline const Tensor4_ddg_Expr<const minus_Tensor4_ddg<A,T,Dim01,Dim23,i,j,k,l>,
  T,Dim01,Dim23,i,j,k,l>
operator-(const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a)
{
  typedef const minus_Tensor4_ddg<A,T,Dim01,Dim23,i,j,k,l> TensorExpr;
  return Tensor4_ddg_Expr<TensorExpr,T,Dim01,Dim23,i,j,k,l>(TensorExpr(a));
}
