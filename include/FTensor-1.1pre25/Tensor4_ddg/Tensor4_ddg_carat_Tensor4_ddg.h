/* This file has all of the declarations for expressions like
   Tensor4_ddg^Tensor4_ddg, yielding a Tensor4_ddg. */

/* A(i,j,k,l)*B(j,l,m,n) */

template<class A, class B, class T, class U, int Dim, int Dim23,
  char i, char j, char k, char l, char m, char n>
class Tensor4_ddg_carat_Tensor4_ddg_13
{
  const Tensor4_ddg_Expr<A,T,Dim,Dim,i,j,k,l> iterA;
  const Tensor4_ddg_Expr<B,U,Dim,Dim23,j,l,m,n> iterB;

  template<int Current_Dim0, int Current_Dim1>
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const int N3, const int N4,
		       const Number<Current_Dim0> &ND0,
		       const Number<Current_Dim1> &ND1)  const
  {
    return iterA(Current_Dim0-1,N1,Current_Dim1-1,N2)
      *iterB(Current_Dim0-1,Current_Dim1-1,N3,N4)
      + eval(N1,N2,N3,N4,Number<Current_Dim0-1>(),Number<Current_Dim1>());
  }
  template<int Current_Dim1>
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const int N3, const int N4,
		       const Number<1> &ND0,
		       const Number<Current_Dim1> &ND1)  const
  {
    return iterA(0,N1,Current_Dim1-1,N2)*iterB(0,Current_Dim1-1,N3,N4)
      + eval(N1,N2,N3,N4,Number<Dim>(),Number<Current_Dim1-1>());
  }
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const int N3, const int N4,
		       const Number<1> &ND0,
		       const Number<1> &ND1)  const
  {
    return iterA(0,N1,0,N2)*iterB(0,0,N3,N4);
  }
public:
  Tensor4_ddg_carat_Tensor4_ddg_13
  (const Tensor4_ddg_Expr<A,T,Dim,Dim,i,j,k,l> &a,
   const Tensor4_ddg_Expr<B,U,Dim,Dim23,j,l,m,n> &b): iterA(a), iterB(b) {}
  typename promote<T,U>::V operator()(const int N1, const int N2,
			     const int N3, const int N4) const
  {
    return eval(N1,N2,N3,N4,Number<Dim>(),Number<Dim>());
  }
};

template<class A, class B, class T, class U, int Dim, int Dim23,
  char i, char j, char k, char l, char m, char n>
inline const Tensor4_ddg_Expr
<const Tensor4_ddg_carat_Tensor4_ddg_13<A,B,T,U,Dim,Dim23,i,j,k,l,m,n>,
  typename promote<T,U>::V,Dim,Dim23,i,k,m,n>
operator^(const Tensor4_ddg_Expr<A,T,Dim,Dim,i,j,k,l> &a,
	  const Tensor4_ddg_Expr<B,U,Dim,Dim23,j,l,m,n> &b)
{
  typedef const Tensor4_ddg_carat_Tensor4_ddg_13<A,B,T,U,Dim,Dim23,i,j,k,l,m,n>
    TensorExpr;
  return Tensor4_ddg_Expr<TensorExpr,typename promote<T,U>::V,Dim,Dim23,i,k,m,n>
    (TensorExpr(a,b));
}
