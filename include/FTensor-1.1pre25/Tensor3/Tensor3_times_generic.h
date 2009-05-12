/* Multiplies a Tensor3 with a generic, yielding a Tensor3. */

/* A(i,j,k)*generic */

template<class A, class T, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k>
class Tensor3_times_generic
{
  const Tensor3_Expr<A,T,Dim0,Dim1,Dim2,i,j,k> iterA;
  const U d;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2, const int N3) const
  {
    return iterA(N1,N2,N3)*d;
  }

  Tensor3_times_generic(const Tensor3_Expr<A,T,Dim0,Dim1,Dim2,i,j,k> &a,
			const U &d0): iterA(a), d(d0) {}
};

template<class A, class T, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k>
inline const Tensor3_Expr<const Tensor3_times_generic<A,T,U,Dim0,Dim1,Dim2,i,j,k>,
  typename promote<T,U>::V,Dim0,Dim1,Dim2,i,j,k>
operator*(const Tensor3_Expr<A,T,Dim0,Dim1,Dim2,i,j,k> &a, const U &d0)
{
  typedef const Tensor3_times_generic<A,T,U,Dim0,Dim1,Dim2,i,j,k>
    TensorExpr;
  return Tensor3_Expr<TensorExpr,typename promote<T,U>::V,Dim0,Dim1,Dim2,i,j,k>
    (TensorExpr(a,d0));
}

/* generic*A(i,j,k) */

template<class A, class T, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k>
inline const Tensor3_Expr<const Tensor3_times_generic<A,T,U,Dim0,Dim1,Dim2,i,j,k>,
  typename promote<T,U>::V,Dim0,Dim1,Dim2,i,j,k>
operator*(const U &d0, const Tensor3_Expr<A,T,Dim0,Dim1,Dim2,i,j,k> &a)
{
  typedef const Tensor3_times_generic<A,T,U,Dim0,Dim1,Dim2,i,j,k>
    TensorExpr;
  return Tensor3_Expr<TensorExpr,typename promote<T,U>::V,Dim0,Dim1,Dim2,i,j,k>
    (TensorExpr(a,d0));
}


