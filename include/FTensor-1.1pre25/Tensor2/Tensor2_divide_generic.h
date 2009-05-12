/* Divides a Tensor2 by a generic, yielding a Tensor2. */

template<class A, class T, class U, int Dim0, int Dim1, char i, char j>
class Tensor2_divide_generic
{
  const Tensor2_Expr<A,T,Dim0,Dim1,i,j> iterA;
  const U d;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return iterA(N1,N2)/d;
  }

  Tensor2_divide_generic(const Tensor2_Expr<A,T,Dim0,Dim1,i,j> &a,
			 const U &d0):
    iterA(a), d(d0) {}
};

template<class A, class T, class U, int Dim0, int Dim1, char i, char j>
inline const Tensor2_Expr<const Tensor2_divide_generic<A,T,U,Dim0,Dim1,i,j>,
  typename promote<T,U>::V,Dim0,Dim1,i,j>
operator/(const Tensor2_Expr<A,T,Dim0,Dim1,i,j> &a, const U &d0)
{
  typedef const Tensor2_divide_generic<A,T,U,Dim0,Dim1,i,j> TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim0,Dim1,i,j>
    (TensorExpr(a,d0));
}

