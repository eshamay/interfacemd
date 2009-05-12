/* This file has all of the declarations for expressions like
   Tensor4_ddg*Tensor1 and Tensor1*Tensor4_ddg, yielding a
   Tensor3_dg. */

/* A(i,j,k,l)*B(k)->Tensor3_dg */

template<class A, class B, class T, class U, int Dim01, int Dim23,
  char i, char j, char k, char l>
class Tensor4_ddg_times_Tensor1_2
{
  const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> iterA;
  const Tensor1_Expr<B,U,Dim23,k> iterB;

  template<int Current_Dim>
  typename promote<T,U>::V eval(const int N1, const int N2, const int N3,
		       const Number<Current_Dim> &N)  const
  {
    return iterA(N1,N2,Current_Dim-1,N3)*iterB(Current_Dim-1)
      + eval(N1,N2,N3,Number<Current_Dim-1>());
  }
  typename promote<T,U>::V eval(const int N1, const int N2, const int N3,
		       const Number<1> &N) const
  {
    return iterA(N1,N2,0,N3)*iterB(0);
  }
public:
  Tensor4_ddg_times_Tensor1_2
  (const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a,
   const Tensor1_Expr<B,U,Dim23,k> &b): iterA(a), iterB(b) {}
  typename promote<T,U>::V operator()(const int N1, const int N2, const int N3) const
  {
    return eval(N1,N2,N3,Number<Dim23>());
  }
};

template<class A, class B, class T, class U, int Dim01, int Dim23,
  char i, char j, char k, char l>
inline const Tensor3_dg_Expr
<const Tensor4_ddg_times_Tensor1_2<A,B,T,U,Dim01,Dim23,i,j,k,l>,
  typename promote<T,U>::V,Dim01,Dim23,i,j,l>
operator*(const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a,
	  const Tensor1_Expr<B,U,Dim23,k> &b)
{
  typedef const Tensor4_ddg_times_Tensor1_2<A,B,T,U,Dim01,Dim23,i,j,k,l>
    TensorExpr;
  return Tensor3_dg_Expr<TensorExpr,typename promote<T,U>::V,Dim01,Dim23,i,j,l>
    (TensorExpr(a,b));
}

/* B(k)*A(i,j,k,l)->Tensor3_dg */

template<class A, class B, class T, class U, int Dim01, int Dim23,
  char i, char j, char k, char l>
inline const Tensor3_dg_Expr
<const Tensor4_ddg_times_Tensor1_2<A,B,T,U,Dim01,Dim23,i,j,k,l>,
  typename promote<T,U>::V,Dim01,Dim23,i,j,l>
operator*(const Tensor1_Expr<B,U,Dim23,k> &b,
	  const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a)
{
  typedef const Tensor4_ddg_times_Tensor1_2<A,B,T,U,Dim01,Dim23,i,j,k,l>
    TensorExpr;
  return Tensor3_dg_Expr<TensorExpr,typename promote<T,U>::V,Dim01,Dim23,i,j,l>
    (TensorExpr(a,b));
}

