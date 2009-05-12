/* Multiply a Tensor1 and a Tensor3_dg together but don't contract, yielding a
   Tensor3_dg. */

/* A(i,j,k) & B(k) -> Tensor3_dg */

template<class A, class B, class T, class U, int Dim01, int Dim2,
  char i, char j, char k>
class Tensor3_dg_and_Tensor1
{
  const Tensor3_dg_Expr<A,T,Dim01,Dim2,i,j,k> iterA;
  const Tensor1_Expr<B,U,Dim2,k> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2, const int N3) const
  {
    return iterA(N1,N2,N3)*iterB(N3);
  }

  Tensor3_dg_and_Tensor1(const Tensor3_dg_Expr<A,T,Dim01,Dim2,i,j,k> &a,
			 const Tensor1_Expr<B,U,Dim2,k> &b):
    iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim01, int Dim2,
  char i, char j, char k>
inline const Tensor3_dg_Expr<const Tensor3_dg_and_Tensor1<A,B,T,U,Dim01,Dim2,i,j,k>,
  typename promote<T,U>::V,Dim01,Dim2,i,j,k>
operator&(const Tensor3_dg_Expr<A,T,Dim01,Dim2,i,j,k> &a,
	   const Tensor1_Expr<B,U,Dim2,k> &b)
{
  typedef const Tensor3_dg_and_Tensor1<A,B,T,U,Dim01,Dim2,i,j,k> TensorExpr;
  return Tensor3_dg_Expr<TensorExpr,typename promote<T,U>::V,Dim01,Dim2,i,j,k>
    (TensorExpr(a,b));
}

/* B(k) & A(i,j,k) -> Tensor3_dg */

template<class A, class B, class T, class U, int Dim01, int Dim2,
  char i, char j, char k>
inline const Tensor3_dg_Expr<const Tensor3_dg_and_Tensor1<A,B,T,U,Dim01,Dim2,i,j,k>,
  typename promote<T,U>::V,Dim01,Dim2,i,j,k>
operator&(const Tensor1_Expr<B,U,Dim2,k> &b,
	   const Tensor3_dg_Expr<A,T,Dim01,Dim2,i,j,k> &a)
{
  typedef const Tensor3_dg_and_Tensor1<A,B,T,U,Dim01,Dim2,i,j,k> TensorExpr;
  return Tensor3_dg_Expr<TensorExpr,typename promote<T,U>::V,Dim01,Dim2,i,j,k>
    (TensorExpr(a,b));
}
