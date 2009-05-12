/* Unary minus operator. */

template<class A, class T, int Dim01, int Dim2, char i, char j, char k>
class minus_Tensor3_dg
{
  const Tensor3_dg_Expr<A,T,Dim01,Dim2,i,j,k> iterA;
public:
  T operator()(const int N1, const int N2, const int N3) const
  {
    return -iterA(N1,N2,N3);
  }

  minus_Tensor3_dg(const Tensor3_dg_Expr<A,T,Dim01,Dim2,i,j,k> &a): iterA(a) {}
};

template<class A, class T, int Dim01, int Dim2, char i, char j, char k>
inline const Tensor3_dg_Expr<const minus_Tensor3_dg<A,T,Dim01,Dim2,i,j,k>,
  T,Dim01,Dim2,i,j,k>
operator-(const Tensor3_dg_Expr<A,T,Dim01,Dim2,i,j,k> &a)
{
    typedef const minus_Tensor3_dg<A,T,Dim01,Dim2,i,j,k> TensorExpr;
    return Tensor3_dg_Expr<TensorExpr,T,Dim01,Dim2,i,j,k>(TensorExpr(a));
}

