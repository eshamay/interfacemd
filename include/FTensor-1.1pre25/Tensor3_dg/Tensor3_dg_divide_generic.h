/* Divides a Tensor3_dg by a generic, yielding a Tensor3_dg. */

template<class A, class T, class U, int Dim01, int Dim2, char i, char j,char k>
class Tensor3_dg_divide_generic
{
  const Tensor3_dg_Expr<A,T,Dim01,Dim2,i,j,k> iterA;
  const U d;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2, const int N3) const
  {
    return iterA(N1,N2,N3)/d;
  }

  Tensor3_dg_divide_generic(const Tensor3_dg_Expr<A,T,Dim01,Dim2,i,j,k> &a,
			   const U &d0): iterA(a), d(d0) {}
};

/* A(i,j,k)/d0->Tensor3_dg */

template<class A, class T, class U, int Dim01, int Dim2, char i, char j,char k>
inline const Tensor3_dg_Expr<const Tensor3_dg_divide_generic<A,T,U,Dim01,Dim2,i,j,k>,
  typename promote<T,U>::V,Dim01,Dim2,i,j,k>
operator/(const Tensor3_dg_Expr<A,T,Dim01,Dim2,i,j,k> &a, const U &d0)
{
  typedef const Tensor3_dg_divide_generic<A,T,U,Dim01,Dim2,i,j,k> TensorExpr;
  return Tensor3_dg_Expr<TensorExpr,typename promote<T,U>::V,Dim01,Dim2,i,j,k>
    (TensorExpr(a,d0));
}
