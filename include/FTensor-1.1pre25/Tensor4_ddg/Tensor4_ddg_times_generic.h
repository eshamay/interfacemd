/* Multiplies a Tensor4_ddg with a generic, yielding a
   Tensor4_ddg. */

template<class A, class T, class U, int Dim01, int Dim23,
  char i, char j, char k, char l>
class Tensor4_ddg_times_generic
{
  const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> iterA;
  const U d;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2, const int N3,
			     const int N4) const
  {
    return iterA(N1,N2,N3,N4)*d;
  }

  Tensor4_ddg_times_generic(const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a,
				  const U &d0): iterA(a), d(d0) {}
};

template<class A, class T, class U, int Dim01, int Dim23,
  char i, char j, char k, char l>
inline const Tensor4_ddg_Expr
<const Tensor4_ddg_times_generic<A,T,U,Dim01,Dim23,i,j,k,l>,typename promote<T,U>::V,Dim01,Dim23,i,j,k,l>
operator*(const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a, const U &d0)
{
  typedef const Tensor4_ddg_times_generic<A,T,U,Dim01,Dim23,i,j,k,l>
    TensorExpr;
  return Tensor4_ddg_Expr<TensorExpr,typename promote<T,U>::V,Dim01,Dim23,i,j,k,l>
    (TensorExpr(a,d0));
}

template<class A, class T, class U, int Dim01, int Dim23,
  char i, char j, char k, char l>
inline const Tensor4_ddg_Expr
<const Tensor4_ddg_times_generic<A,T,U,Dim01,Dim23,i,j,k,l>,typename promote<T,U>::V,Dim01,Dim23,i,j,k,l>
operator*(const U &d0, const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a)
{
  typedef const Tensor4_ddg_times_generic<A,T,U,Dim01,Dim23,i,j,k,l>
    TensorExpr;
  return Tensor4_ddg_Expr<TensorExpr,typename promote<T,U>::V,Dim01,Dim23,i,j,k,l>
    (TensorExpr(a,d0));
}
