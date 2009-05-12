/* Adds two Tensor1's together yielding a Tensor2_symmetric. */

template<class A, class B, class T, class U, int Dim, char i, char j>
class Tensor1_or_Tensor1
{
  const Tensor1_Expr<A,T,Dim,i> iterA;
  const Tensor1_Expr<B,U,Dim,j> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return iterA(N1)+iterB(N2);
  }

  Tensor1_or_Tensor1(const Tensor1_Expr<A,T,Dim,i> &a,
			const Tensor1_Expr<B,U,Dim,j> &b)
    : iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim, char i, char j>
inline const Tensor2_symmetric_Expr<const Tensor1_or_Tensor1<A,B,T,U,Dim,i,j>,
  typename promote<T,U>::V,Dim,i,j>
operator||(const Tensor1_Expr<A,T,Dim,i> &a, const Tensor1_Expr<B,U,Dim,j> &b)
{
  typedef const Tensor1_or_Tensor1<A,B,T,U,Dim,i,j> TensorExpr;
  return Tensor2_symmetric_Expr<TensorExpr,typename promote<T,U>::V,Dim,i,j>
    (TensorExpr(a,b));
}
