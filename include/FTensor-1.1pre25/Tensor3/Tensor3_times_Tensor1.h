/* Declarations for expressions like Tensor3*Tensor1 -> Tensor2 */

/* A(i,j,k)*B(j) -> Tensor2 */

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k>
class Tensor3_times_Tensor1_1
{
  const Tensor3_Expr<A,T,Dim0,Dim1,Dim2,i,j,k> iterA;
  const Tensor1_Expr<B,U,Dim1,j> iterB;

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
  Tensor3_times_Tensor1_1(const Tensor3_Expr<A,T,Dim0,Dim1,Dim2,i,j,k> &a,
			const Tensor1_Expr<B,U,Dim1,j> &b): iterA(a), iterB(b) {}
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return eval(N1,N2,Number<Dim1>());
  }
};

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k>
inline const Tensor2_Expr
<const Tensor3_times_Tensor1_1<A,B,T,U,Dim0,Dim1,Dim2,i,j,k>,
  typename promote<T,U>::V,Dim0,Dim2,i,k>
operator*(const Tensor3_Expr<A,T,Dim0,Dim1,Dim2,i,j,k> &a,
	  const Tensor1_Expr<B,U,Dim1,j> &b)
{
  typedef const Tensor3_times_Tensor1_1<A,B,T,U,Dim0,Dim1,Dim2,i,j,k>
    TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim0,Dim2,i,k>
    (TensorExpr(a,b));
}

/* B(j)*A(i,j,k) -> Tensor2 */

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k>
inline const Tensor2_Expr
<const Tensor3_times_Tensor1_1<A,B,T,U,Dim0,Dim1,Dim2,i,j,k>,
  typename promote<T,U>::V,Dim0,Dim2,i,k>
operator*(const Tensor1_Expr<B,U,Dim1,j> &b,
	  const Tensor3_Expr<A,T,Dim0,Dim1,Dim2,i,j,k> &a)
{
  typedef const Tensor3_times_Tensor1_1<A,B,T,U,Dim0,Dim1,Dim2,i,j,k>
    TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim0,Dim2,i,k>
    (TensorExpr(a,b));
}

/* A(i,j,k)*B(i) -> Tensor2 */

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k>
class Tensor3_times_Tensor1_0
{
  const Tensor3_Expr<A,T,Dim0,Dim1,Dim2,i,j,k> iterA;
  const Tensor1_Expr<B,U,Dim0,i> iterB;

  template<int Current_Dim>
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<Current_Dim> &N)  const
  {
    return iterA(Current_Dim-1,N1,N2)*iterB(Current_Dim-1)
      + eval(N1,N2,Number<Current_Dim-1>());
  }
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<1> &N) const
  {
    return iterA(0,N1,N2)*iterB(0);
  }
public:
  Tensor3_times_Tensor1_0(const Tensor3_Expr<A,T,Dim0,Dim1,Dim2,i,j,k> &a,
			  const Tensor1_Expr<B,U,Dim0,i> &b)
    : iterA(a), iterB(b) {}
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return eval(N1,N2,Number<Dim0>());
  }
};

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k>
inline const Tensor2_Expr
<const Tensor3_times_Tensor1_0<A,B,T,U,Dim0,Dim1,Dim2,i,j,k>,
  typename promote<T,U>::V,Dim1,Dim2,j,k>
operator*(const Tensor3_Expr<A,T,Dim0,Dim1,Dim2,i,j,k> &a,
	  const Tensor1_Expr<B,U,Dim0,i> &b)
{
  typedef const Tensor3_times_Tensor1_0<A,B,T,U,Dim0,Dim1,Dim2,i,j,k>
    TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim1,Dim2,j,k>
    (TensorExpr(a,b));
}

/* B(i)*A(i,j,k) -> Tensor2 */

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k>
inline const Tensor2_Expr
<const Tensor3_times_Tensor1_0<A,B,T,U,Dim0,Dim1,Dim2,i,j,k>,
  typename promote<T,U>::V,Dim1,Dim2,j,k>
operator*(const Tensor1_Expr<B,U,Dim0,i> &b,
	  const Tensor3_Expr<A,T,Dim0,Dim1,Dim2,i,j,k> &a)
{
  typedef const Tensor3_times_Tensor1_0<A,B,T,U,Dim0,Dim1,Dim2,i,j,k>
    TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim1,Dim2,j,k>
    (TensorExpr(a,b));
}

