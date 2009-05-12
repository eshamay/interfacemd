/* Unary minus operator. */

template<class A, class T, int Dim, char i, char j>
class minus_Tensor2_symmetric
{
public:
  const Tensor2_symmetric_Expr<A,T,Dim,i,j> iterA;
public:
  T operator()(const int N1, const int N2) const
  {
    return -iterA(N1,N2);
  }
  minus_Tensor2_symmetric(const Tensor2_symmetric_Expr<A,T,Dim,i,j> &a):
    iterA(a) {}
};

template<class A, class T, int Dim, char i, char j>
inline const Tensor2_symmetric_Expr<const minus_Tensor2_symmetric<A,T,Dim,i,j>,
  T,Dim,i,j>
operator-(const Tensor2_symmetric_Expr<A,T,Dim,i,j> &a)
{
  typedef const minus_Tensor2_symmetric<A,T,Dim,i,j> TensorExpr;
  return Tensor2_symmetric_Expr<TensorExpr,T,Dim,i,j>(TensorExpr(a));
}
