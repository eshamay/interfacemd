/* Adds two Tensor1's together, yielding a Tensor1. */

template<class A, class B, class T, class U, int Dim, char i>
class Tensor1_plus_Tensor1
{
  const Tensor1_Expr<A,T,Dim,i> iterA;
  const Tensor1_Expr<B,U,Dim,i> iterB;
public:
  typename promote<T,U>::V operator()(const int N) const
  {
    return iterA(N)+iterB(N);
  }

  Tensor1_plus_Tensor1(const Tensor1_Expr<A,T,Dim,i> &a,
		       const Tensor1_Expr<B,U,Dim,i> &b):
    iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, char i, int Dim>
inline const Tensor1_Expr<const Tensor1_plus_Tensor1<A,B,T,U,Dim,i>,
  typename promote<T,U>::V,Dim,i>
operator+(const Tensor1_Expr<A,T,Dim,i> &a, const Tensor1_Expr<B,U,Dim,i> &b)
{
  typedef const Tensor1_plus_Tensor1<A,B,T,U,Dim,i> TensorExpr;
  return Tensor1_Expr<TensorExpr,typename promote<T,U>::V,Dim,i>
    (TensorExpr(a,b));
}
