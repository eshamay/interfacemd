/* Divide a Tensor4_ddg by a Tensor2_symmetric without contracting, yielding a
   Tensor4_ddg. */

/* A(i,j,k,l) % B(i,j) -> Tensor4_ddg */

template<class A, class B, class T, class U, int Dim01, int Dim23,
  char i, char j, char k, char l>
class Tensor4_ddg_mod_Tensor2_symmetric_01
{
  const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> iterA;
  const Tensor2_symmetric_Expr<B,U,Dim01,i,j> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2, const int N3,
		    const int N4) const
  {
    return iterA(N1,N2,N3,N4)/iterB(N1,N2);
  }

  Tensor4_ddg_mod_Tensor2_symmetric_01
  (const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a,
   const Tensor2_symmetric_Expr<B,U,Dim01,i,j> &b): iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim01, int Dim23,
  char i, char j, char k, char l>
inline const Tensor4_ddg_Expr
<const Tensor4_ddg_mod_Tensor2_symmetric_01<A,B,T,U,Dim01,Dim23,i,j,k,l>,
  typename promote<T,U>::V,Dim01,Dim23,i,j,k,l>
operator%(const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a,
	  const Tensor2_symmetric_Expr<B,U,Dim01,i,j> &b)
{
  typedef const Tensor4_ddg_mod_Tensor2_symmetric_01<A,B,T,U,Dim01,Dim23,i,j,k,l>
    TensorExpr;
  return Tensor4_ddg_Expr<TensorExpr,typename promote<T,U>::V,Dim01,Dim23,i,j,k,l>
    (TensorExpr(a,b));
}

/* B(i,j) % A(i,j,k,l) -> Tensor4_ddg */

template<class A, class B, class T, class U, int Dim01, int Dim23,
  char i, char j, char k, char l>
inline const Tensor4_ddg_Expr
<const Tensor4_ddg_mod_Tensor2_symmetric_01<A,B,T,U,Dim01,Dim23,i,j,k,l>,
  typename promote<T,U>::V,Dim01,Dim23,i,j,k,l>
operator%(const Tensor2_symmetric_Expr<B,U,Dim01,i,j> &b,
	  const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a)
{
  typedef const Tensor4_ddg_mod_Tensor2_symmetric_01<A,B,T,U,Dim01,Dim23,i,j,k,l>
    TensorExpr;
  return Tensor4_ddg_Expr<TensorExpr,typename promote<T,U>::V,Dim01,Dim23,i,j,k,l>
    (TensorExpr(a,b));
}

/* A(i,j,k,l) % B(k,l) -> Tensor4_ddg */

template<class A, class B, class T, class U, int Dim01, int Dim23,
  char i, char j, char k, char l>
class Tensor4_ddg_mod_Tensor2_symmetric_23
{
  const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> iterA;
  const Tensor2_symmetric_Expr<B,U,Dim01,k,l> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2, const int N3,
		    const int N4) const
  {
    return iterA(N1,N2,N3,N4)/iterB(N3,N4);
  }

  Tensor4_ddg_mod_Tensor2_symmetric_23
  (const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a,
   const Tensor2_symmetric_Expr<B,U,Dim01,k,l> &b): iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim01, int Dim23,
  char i, char j, char k, char l>
inline const Tensor4_ddg_Expr
<const Tensor4_ddg_mod_Tensor2_symmetric_23<A,B,T,U,Dim01,Dim23,i,j,k,l>,
  typename promote<T,U>::V,Dim01,Dim23,i,j,k,l>
operator%(const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a,
	  const Tensor2_symmetric_Expr<B,U,Dim01,k,l> &b)
{
  typedef const Tensor4_ddg_mod_Tensor2_symmetric_23<A,B,T,U,Dim01,Dim23,i,j,k,l>
    TensorExpr;
  return Tensor4_ddg_Expr<TensorExpr,typename promote<T,U>::V,Dim01,Dim23,i,j,k,l>
    (TensorExpr(a,b));
}

/* B(k,l) % A(i,j,k,l) -> Tensor4_ddg */

template<class A, class B, class T, class U, int Dim01, int Dim23,
  char i, char j, char k, char l>
inline const Tensor4_ddg_Expr
<const Tensor4_ddg_mod_Tensor2_symmetric_23<A,B,T,U,Dim01,Dim23,i,j,k,l>,
  typename promote<T,U>::V,Dim01,Dim23,i,j,k,l>
operator%(const Tensor2_symmetric_Expr<B,U,Dim01,k,l> &b,
	  const Tensor4_ddg_Expr<A,T,Dim01,Dim23,i,j,k,l> &a)
{
  typedef const Tensor4_ddg_mod_Tensor2_symmetric_23<A,B,T,U,Dim01,Dim23,i,j,k,l>
    TensorExpr;
  return Tensor4_ddg_Expr<TensorExpr,typename promote<T,U>::V,Dim01,Dim23,i,j,k,l>
    (TensorExpr(a,b));
}
