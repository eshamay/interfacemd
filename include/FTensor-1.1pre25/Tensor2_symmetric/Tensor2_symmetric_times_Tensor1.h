/* This file has all of the declarations for expressions like
   Tensor2_symmetric*Tensor1 and Tensor1*Tensor2_symmetric, yielding a
   Tensor1 or Tensor3_dg. */

/* A(i,j)*B(j) -> Tensor1 */

template<class A, class B, class T, class U, int Dim, char i,char j>
class Tensor2_symmetric_times_Tensor1_1
{
  const Tensor2_symmetric_Expr<A,T,Dim,i,j> iterA;
  const Tensor1_Expr<B,U,Dim,j> iterB;

  template<int Current_Dim>
  typename promote<T,U>::V eval(const int N1, const Number<Current_Dim> &N) const
  {
    return iterA(N1,Current_Dim-1)*iterB(Current_Dim-1)
      + eval(N1,Number<Current_Dim-1>());
  }
  typename promote<T,U>::V eval(const int N1, const Number<1> &N) const
  {
    return iterA(N1,0)*iterB(0);
  }
public:
  Tensor2_symmetric_times_Tensor1_1
  (const Tensor2_symmetric_Expr<A,T,Dim,i,j> &a,
   const Tensor1_Expr<B,U,Dim,j> &b): iterA(a), iterB(b) {}
  typename promote<T,U>::V operator()(const int N1) const
  {
    return eval(N1,Number<Dim>());
  }
};

template<class A, class B, class T, class U, int Dim, char i,char j>
inline const Tensor1_Expr<const Tensor2_symmetric_times_Tensor1_1<A,B,T,U,Dim,i,j>,
  typename promote<T,U>::V,Dim,i>
operator*(const Tensor2_symmetric_Expr<A,T,Dim,i,j> &a,
	  const Tensor1_Expr<B,U,Dim,j> &b)
{
  typedef const Tensor2_symmetric_times_Tensor1_1<A,B,T,U,Dim,i,j> TensorExpr;
  return Tensor1_Expr<TensorExpr,typename promote<T,U>::V,Dim,i>(TensorExpr(a,b));
}

/* A(j,i)*B(j) -> Tensor1 */

template<class A, class B, class T, class U, int Dim, char i,char j>
class Tensor2_symmetric_times_Tensor1_0
{
  const Tensor2_symmetric_Expr<A,T,Dim,j,i> iterA;
  const Tensor1_Expr<B,U,Dim,j> iterB;

  template<int Current_Dim>
  typename promote<T,U>::V eval(const int N1, const Number<Current_Dim> &N) const
  {
    return iterA(Current_Dim-1,N1)*iterB(Current_Dim-1)
      + eval(N1,Number<Current_Dim-1>());
  }
  typename promote<T,U>::V eval(const int N1, const Number<1> &N) const
  {
    return iterA(0,N1)*iterB(0);
  }
public:
  Tensor2_symmetric_times_Tensor1_0
  (const Tensor2_symmetric_Expr<A,T,Dim,j,i> &a,
   const Tensor1_Expr<B,U,Dim,j> &b): iterA(a), iterB(b) {}
  typename promote<T,U>::V operator()(const int N1) const
  {
    return eval(N1,Number<Dim>());
  }
};

template<class A, class B, class T, class U, int Dim, char i,char j>
inline const Tensor1_Expr<const Tensor2_symmetric_times_Tensor1_0<A,B,T,U,Dim,i,j>,
  typename promote<T,U>::V,Dim,i>
operator*(const Tensor2_symmetric_Expr<A,T,Dim,j,i> &a,
	  const Tensor1_Expr<B,U,Dim,j> &b)
{
  typedef const Tensor2_symmetric_times_Tensor1_0<A,B,T,U,Dim,i,j> TensorExpr;
  return Tensor1_Expr<TensorExpr,typename promote<T,U>::V,Dim,i>(TensorExpr(a,b));
}

/* B(j)*A(i,j) -> Tensor1 */

template<class A, class B, class T, class U, int Dim, char i,char j>
inline const Tensor1_Expr<const Tensor2_symmetric_times_Tensor1_1<A,B,T,U,Dim,i,j>,
  typename promote<T,U>::V,Dim,i>
operator*(const Tensor1_Expr<B,U,Dim,j> &b,
	  const Tensor2_symmetric_Expr<A,T,Dim,i,j> &a)
{
  typedef const Tensor2_symmetric_times_Tensor1_1<A,B,T,U,Dim,i,j> TensorExpr;
  return Tensor1_Expr<TensorExpr,typename promote<T,U>::V,Dim,i>(TensorExpr(a,b));
}

/* B(j)*A(j,i) -> Tensor1 */

template<class A, class B, class T, class U, int Dim, char i,char j>
inline const Tensor1_Expr<const Tensor2_symmetric_times_Tensor1_0<A,B,T,U,Dim,i,j>,
  typename promote<T,U>::V,Dim,i>
operator*(const Tensor1_Expr<B,U,Dim,j> &b,
	  const Tensor2_symmetric_Expr<A,T,Dim,j,i> &a)
{
  typedef const Tensor2_symmetric_times_Tensor1_0<A,B,T,U,Dim,i,j> TensorExpr;
  return Tensor1_Expr<TensorExpr,typename promote<T,U>::V,Dim,i>(TensorExpr(a,b));
}

/* A(i,j)*B(k) -> Tensor3_dg */

template<class A, class B, class T, class U, int Dim, int Dim2,
  char i, char j, char k>
class Tensor2_symmetric_times_Tensor1
{
  const Tensor2_symmetric_Expr<A,T,Dim,i,j> iterA;
  const Tensor1_Expr<B,U,Dim2,k> iterB;
public:
  Tensor2_symmetric_times_Tensor1(const Tensor2_symmetric_Expr<A,T,Dim,i,j> &a,
			const Tensor1_Expr<B,U,Dim2,k> &b):
    iterA(a), iterB(b) {}
  typename promote<T,U>::V operator()(const int N1, const int N2, const int N3) const
  {
    return iterA(N1,N2)*iterB(N3);
  }
};

template<class A, class B, class T, class U, int Dim, int Dim2,
  char i, char j, char k>
inline const Tensor3_dg_Expr
<const Tensor2_symmetric_times_Tensor1<A,B,T,U,Dim,Dim2,i,j,k>,
  typename promote<T,U>::V,Dim,Dim2,i,j,k>
operator*(const Tensor2_symmetric_Expr<A,T,Dim,i,j> &a,
	  const Tensor1_Expr<B,U,Dim2,k> &b)
{
  typedef const Tensor2_symmetric_times_Tensor1<A,B,T,U,Dim,Dim2,i,j,k>
    TensorExpr;
  return Tensor3_dg_Expr<TensorExpr,typename promote<T,U>::V,Dim,Dim2,i,j,k>
    (TensorExpr(a,b));
}

/* B(k)*A(i,j) -> Tensor3_dg */

template<class A, class B, class T, class U, int Dim, int Dim2,
  char i, char j, char k>
inline const Tensor3_dg_Expr
<const Tensor2_symmetric_times_Tensor1<A,B,T,U,Dim,Dim2,i,j,k>,
  typename promote<T,U>::V,Dim,Dim2,i,j,k>
operator*(const Tensor1_Expr<B,U,Dim2,k> &b,
	  const Tensor2_symmetric_Expr<A,T,Dim,i,j> &a)
{
  typedef const Tensor2_symmetric_times_Tensor1<A,B,T,U,Dim,Dim2,i,j,k>
    TensorExpr;
  return Tensor3_dg_Expr<TensorExpr,typename promote<T,U>::V,Dim,Dim2,i,j,k>
    (TensorExpr(a,b));
}
