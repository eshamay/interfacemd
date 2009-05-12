/* Unary minus operator. */

template<class A, class T, int Dim0, int Dim1, char i, char j>
class minus_Tensor2
{
  const Tensor2_Expr<A,T,Dim0,Dim1,i,j> iterA;
public:
  T operator()(const int N1, const int N2) const
  {
    return -iterA(N1,N2);
  }

  minus_Tensor2(const Tensor2_Expr<A,T,Dim0,Dim1,i,j> &a): iterA(a) {}
};

template<class A, class T, int Dim0, int Dim1, char i, char j>
inline const Tensor2_Expr<const minus_Tensor2<A,T,Dim0,Dim1,i,j>,T,Dim0,Dim1,i,j>
operator-(const Tensor2_Expr<A,T,Dim0,Dim1,i,j> &a)
{
  typedef const minus_Tensor2<A,T,Dim0,Dim1,i,j> TensorExpr;
  return Tensor2_Expr<TensorExpr,T,Dim0,Dim1,i,j>(TensorExpr(a));
}

