/* Applies an arbitrary function to the Tensor2_Expr. */

template<class A, class B, class T, int Dim0, int Dim1, char i, char j>
class transform_Tensor2
{
  const Tensor2_Expr<A,T,Dim0,Dim1,i,j> iterA;
  B function;
public:
  T operator()(const int N1, const int N2) const
  {
    return function(iterA(N1,N2));
  }

  transform_Tensor2(const Tensor2_Expr<A,T,Dim0,Dim1,i,j> &a, B func):
    iterA(a), function(func) {}
};

template<class A, class B, class T, int Dim0, int Dim1, char i, char j>
inline const Tensor2_Expr<const transform_Tensor2<A,B,T,Dim0,Dim1,i,j>,T,Dim0,Dim1,i,j>
transform(const Tensor2_Expr<A,T,Dim0,Dim1,i,j> &a, B function)
{
  typedef const transform_Tensor2<A,B,T,Dim0,Dim1,i,j> TensorExpr;
  return Tensor2_Expr<TensorExpr,T,Dim0,Dim1,i,j>(TensorExpr(a,function));
}

