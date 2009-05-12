/* Multiply a Tensor2_symmetric and a Tensor3_dg together but don't
   contract, yielding a Tensor3_dg. */

/* A(i,j,k) & B(i,j) -> Tensor3_dg */

template<class A, class B, class T, class U, int Dim01, int Dim2,
  char i, char j, char k>
class Tensor3_dg_and_Tensor2_symmetric
{
  const Tensor3_dg_Expr<A,T,Dim01,Dim2,i,j,k> iterA;
  const Tensor2_symmetric_Expr<B,U,Dim01,i,j> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2, const int N3) const
  {
    return iterA(N1,N2,N3)*iterB(N1,N2);
  }

  Tensor3_dg_and_Tensor2_symmetric
  (const Tensor3_dg_Expr<A,T,Dim01,Dim2,i,j,k> &a,
   const Tensor2_symmetric_Expr<B,U,Dim01,i,j> &b): iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim01, int Dim2,
  char i, char j, char k>
inline const Tensor3_dg_Expr
<const Tensor3_dg_and_Tensor2_symmetric<A,B,T,U,Dim01,Dim2,i,j,k>,
  typename promote<T,U>::V,Dim01,Dim2,i,j,k>
operator&(const Tensor3_dg_Expr<A,T,Dim01,Dim2,i,j,k> &a,
	  const Tensor2_symmetric_Expr<B,U,Dim01,i,j> &b)
{
  typedef const Tensor3_dg_and_Tensor2_symmetric<A,B,T,U,Dim01,Dim2,i,j,k>
    TensorExpr;
  return Tensor3_dg_Expr<TensorExpr,typename promote<T,U>::V,Dim01,Dim2,i,j,k>
    (TensorExpr(a,b));
}

/* B(i,j) & A(i,j,k) -> Tensor3_dg */

template<class A, class B, class T, class U, int Dim01, int Dim2,
  char i, char j, char k>
inline const Tensor3_dg_Expr
<const Tensor3_dg_and_Tensor2_symmetric<A,B,T,U,Dim01,Dim2,i,j,k>,
  typename promote<T,U>::V,Dim01,Dim2,i,j,k>
operator&(const Tensor2_symmetric_Expr<B,U,Dim01,i,j> &b,
	  const Tensor3_dg_Expr<A,T,Dim01,Dim2,i,j,k> &a)
{
  typedef const Tensor3_dg_and_Tensor2_symmetric<A,B,T,U,Dim01,Dim2,i,j,k>
    TensorExpr;
  return Tensor3_dg_Expr<TensorExpr,typename promote<T,U>::V,Dim01,Dim2,i,j,k>
    (TensorExpr(a,b));
}
