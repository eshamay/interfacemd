/* Declares a wrapper class for the unary minus (-) operator. */

template<class A, class T, int Dim, char i>
class minus_Tensor1
{
  const Tensor1_Expr<A,T,Dim,i> iterA;
public:
  T operator()(const int N) const
  {
    return -iterA(N);
  }

  minus_Tensor1(const Tensor1_Expr<A,T,Dim,i> &a): iterA(a) {}
};

template<class A, class T, int Dim, char i>
inline const Tensor1_Expr<const minus_Tensor1<A,T,Dim,i>,T,Dim,i>
operator-(const Tensor1_Expr<A,T,Dim,i> &a)
{
    typedef const minus_Tensor1<A,T,Dim,i> TensorExpr;
    return Tensor1_Expr<TensorExpr,T,Dim,i>(TensorExpr(a));
}
