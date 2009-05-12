/* This file has all of the declarations for expressions like
   Tensor3_dg*Tensor1 and Tensor1*Tensor3_dg, yielding a
   Tensor2_symmetric or Tensor2. */

/* A(i,j,k)*B(k)->Tensor2_symmetric */

template<class A, class B, class T, class U, int Dim01, int Dim2,
  char i, char j, char k>
class Tensor3_dg_times_Tensor1_2
{
  const Tensor3_dg_Expr<A,T,Dim01,Dim2,i,j,k> iterA;
  const Tensor1_Expr<B,U,Dim2,k> iterB;

  template<int Current_Dim>
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<Current_Dim> &N)  const
  {
    return iterA(N1,N2,Current_Dim-1)*iterB(Current_Dim-1)
      + eval(N1,N2,Number<Current_Dim-1>());
  }
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<1> &N) const
  {
    return iterA(N1,N2,0)*iterB(0);
  }
public:
  Tensor3_dg_times_Tensor1_2(const Tensor3_dg_Expr<A,T,Dim01,Dim2,i,j,k> &a,
			     const Tensor1_Expr<B,U,Dim2,k> &b)
    : iterA(a), iterB(b) {}
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return eval(N1,N2,Number<Dim2>());
  }
};

template<class A, class B, class T, class U, int Dim01, int Dim2,
  char i, char j, char k>
inline const Tensor2_symmetric_Expr
<const Tensor3_dg_times_Tensor1_2<A,B,T,U,Dim01,Dim2,i,j,k>,
  typename promote<T,U>::V,Dim01,i,j>
operator*(const Tensor3_dg_Expr<A,T,Dim01,Dim2,i,j,k> &a,
	  const Tensor1_Expr<B,U,Dim2,k> &b)
{
  typedef const Tensor3_dg_times_Tensor1_2<A,B,T,U,Dim01,Dim2,i,j,k>
    TensorExpr;
  return Tensor2_symmetric_Expr<TensorExpr,typename promote<T,U>::V,Dim01,i,j>
    (TensorExpr(a,b));
}

/* B(k)*A(i,j,k)->Tensor2_symmetric */

template<class A, class B, class T, class U, int Dim01, int Dim2,
  char i, char j, char k>
inline const Tensor2_symmetric_Expr
<const Tensor3_dg_times_Tensor1_2<A,B,T,U,Dim01,Dim2,i,j,k>,
  typename promote<T,U>::V,Dim01,i,j>
operator*(const Tensor1_Expr<B,U,Dim2,k> &b,
	  const Tensor3_dg_Expr<A,T,Dim01,Dim2,i,j,k> &a)
{
  typedef const Tensor3_dg_times_Tensor1_2<A,B,T,U,Dim01,Dim2,i,j,k>
    TensorExpr;
  return Tensor2_symmetric_Expr<TensorExpr,typename promote<T,U>::V,Dim01,i,j>
    (TensorExpr(a,b));
}

/* A(i,k,j)*B(k)->Tensor2 */

template<class A, class B, class T, class U, int Dim01, int Dim2,
  char i, char j, char k>
class Tensor3_dg_times_Tensor1_1
{
  const Tensor3_dg_Expr<A,T,Dim01,Dim2,i,k,j> iterA;
  const Tensor1_Expr<B,U,Dim01,k> iterB;

  template<int Current_Dim>
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<Current_Dim> &N)  const
  {
    return iterA(N1,Current_Dim-1,N2)*iterB(Current_Dim-1)
      + eval(N1,N2,Number<Current_Dim-1>());
  }
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<1> &N) const
  {
    return iterA(N1,0,N2)*iterB(0);
  }
public:
  Tensor3_dg_times_Tensor1_1(const Tensor3_dg_Expr<A,T,Dim01,Dim2,i,k,j> &a,
			     const Tensor1_Expr<B,U,Dim01,k> &b)
    : iterA(a), iterB(b) {}
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return eval(N1,N2,Number<Dim01>());
  }
};

template<class A, class B, class T, class U, int Dim01, int Dim2,
  char i, char j, char k>
inline const Tensor2_Expr<const Tensor3_dg_times_Tensor1_1<A,B,T,U,Dim01,Dim2,i,j,k>,
  typename promote<T,U>::V,Dim01,Dim2,i,j>
operator*(const Tensor3_dg_Expr<A,T,Dim01,Dim2,i,k,j> &a,
	  const Tensor1_Expr<B,U,Dim01,k> &b)
{
  typedef const Tensor3_dg_times_Tensor1_1<A,B,T,U,Dim01,Dim2,i,j,k>
    TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim01,Dim2,i,j>
    (TensorExpr(a,b));
}

/* B(k)*A(i,k,j)->Tensor2 */

template<class A, class B, class T, class U, int Dim01, int Dim2,
  char i, char j, char k>
inline const Tensor2_Expr<const Tensor3_dg_times_Tensor1_1<A,B,T,U,Dim01,Dim2,i,j,k>,
  typename promote<T,U>::V,Dim01,Dim2,i,j>
operator*(const Tensor1_Expr<B,U,Dim01,k> &b,
	  const Tensor3_dg_Expr<A,T,Dim01,Dim2,i,k,j> &a)
{
  typedef const Tensor3_dg_times_Tensor1_1<A,B,T,U,Dim01,Dim2,i,j,k>
    TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim01,Dim2,i,j>
    (TensorExpr(a,b));
}

/* A(k,i,j)*B(k)->Tensor2 */

template<class A, class B, class T, class U, int Dim01, int Dim2,
  char i, char j, char k>
class Tensor3_dg_times_Tensor1_0
{
  const Tensor3_dg_Expr<A,T,Dim01,Dim2,k,i,j> iterA;
  const Tensor1_Expr<B,U,Dim01,k> iterB;

  template<int Current_Dim>
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<Current_Dim> &N)  const
  {
    return iterA(N1,Current_Dim-1,N2)*iterB(Current_Dim-1)
      + eval(N1,N2,Number<Current_Dim-1>());
  }
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<1> &N) const
  {
    return iterA(N1,0,N2)*iterB(0);
  }
public:
  Tensor3_dg_times_Tensor1_0(const Tensor3_dg_Expr<A,T,Dim01,Dim2,k,i,j> &a,
			     const Tensor1_Expr<B,U,Dim01,k> &b)
    : iterA(a), iterB(b) {}
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return eval(N1,N2,Number<Dim01>());
  }
};

template<class A, class B, class T, class U, int Dim01, int Dim2,
  char i, char j, char k>
inline const Tensor2_Expr<const Tensor3_dg_times_Tensor1_0<A,B,T,U,Dim01,Dim2,i,j,k>,
  typename promote<T,U>::V,Dim01,Dim2,i,j>
operator*(const Tensor3_dg_Expr<A,T,Dim01,Dim2,k,i,j> &a,
	  const Tensor1_Expr<B,U,Dim01,k> &b)
{
  typedef const Tensor3_dg_times_Tensor1_0<A,B,T,U,Dim01,Dim2,i,j,k>
    TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim01,Dim2,i,j>
    (TensorExpr(a,b));
}

/* B(k)*A(k,i,j)->Tensor2 */

template<class A, class B, class T, class U, int Dim01, int Dim2,
  char i, char j, char k>
inline const Tensor2_Expr<const Tensor3_dg_times_Tensor1_0<A,B,T,U,Dim01,Dim2,i,j,k>,
  typename promote<T,U>::V,Dim01,Dim2,i,j>
operator*(const Tensor1_Expr<B,U,Dim01,k> &b,
	  const Tensor3_dg_Expr<A,T,Dim01,Dim2,k,i,j> &a)
{
  typedef const Tensor3_dg_times_Tensor1_0<A,B,T,U,Dim01,Dim2,i,j,k>
    TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim01,Dim2,i,j>
    (TensorExpr(a,b));
}
