/* This file has all of the declarations for expressions like
   Tensor2*Tensor1 and Tensor1*Tensor2, yielding a Tensor1 or
   Tensor3. */

/* A(i,j)*B(j) -> Tensor1 */

template<class A, class B, class T, class U, int Dim0, int Dim1, char i,char j>
class Tensor2_times_Tensor1_1
{
  const Tensor2_Expr<A,T,Dim0,Dim1,i,j> iterA;
  const Tensor1_Expr<B,U,Dim1,j> iterB;

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
  Tensor2_times_Tensor1_1(const Tensor2_Expr<A,T,Dim0,Dim1,i,j> &a,
			  const Tensor1_Expr<B,U,Dim1,j> &b):
    iterA(a), iterB(b) {}
  typename promote<T,U>::V operator()(const int N1) const
  {
    return eval(N1,Number<Dim1>());
  }
};

template<class A, class B, class T, class U, int Dim0, int Dim1, char i,char j>
inline const Tensor1_Expr<const Tensor2_times_Tensor1_1<A,B,T,U,Dim0,Dim1,i,j>,
  typename promote<T,U>::V,Dim0,i>
operator*(const Tensor2_Expr<A,T,Dim0,Dim1,i,j> &a,
	  const Tensor1_Expr<B,U,Dim1,j> &b)
{
  typedef const Tensor2_times_Tensor1_1<A,B,T,U,Dim0,Dim1,i,j> TensorExpr;
  return Tensor1_Expr<TensorExpr,typename promote<T,U>::V,Dim0,i>(TensorExpr(a,b));
}

/* A(j,i)*B(j) -> Tensor1 */

template<class A, class B, class T, class U, int Dim0, int Dim1, char i,char j>
class Tensor2_times_Tensor1_0
{
  const Tensor2_Expr<A,T,Dim1,Dim0,j,i> iterA;
  const Tensor1_Expr<B,U,Dim1,j> iterB;

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
  Tensor2_times_Tensor1_0(const Tensor2_Expr<A,T,Dim1,Dim0,j,i> &a,
			  const Tensor1_Expr<B,U,Dim1,j> &b):
    iterA(a), iterB(b) {}
  typename promote<T,U>::V operator()(const int N1) const
  {
    return eval(N1,Number<Dim1>());
  }
};

template<class A, class B, class T, class U, int Dim0, int Dim1, char i,char j>
inline const Tensor1_Expr<const Tensor2_times_Tensor1_0<A,B,T,U,Dim0,Dim1,i,j>,
  typename promote<T,U>::V,Dim0,i>
operator*(const Tensor2_Expr<A,T,Dim1,Dim0,j,i> &a,
	  const Tensor1_Expr<B,U,Dim1,j> &b)
{
  typedef const Tensor2_times_Tensor1_0<A,B,T,U,Dim0,Dim1,i,j> TensorExpr;
  return Tensor1_Expr<TensorExpr,typename promote<T,U>::V,Dim0,i>(TensorExpr(a,b));
}

/* B(j)*A(i,j) -> Tensor1 */

template<class A, class B, class T, class U, int Dim0, int Dim1, char i,char j>
inline const Tensor1_Expr<const Tensor2_times_Tensor1_1<A,B,T,U,Dim0,Dim1,i,j>,
  typename promote<T,U>::V,Dim0,i>
operator*(const Tensor1_Expr<B,U,Dim1,j> &b,
	  const Tensor2_Expr<A,T,Dim0,Dim1,i,j> &a)
{
  typedef const Tensor2_times_Tensor1_1<A,B,T,U,Dim0,Dim1,i,j> TensorExpr;
  return Tensor1_Expr<TensorExpr,typename promote<T,U>::V,Dim0,i>(TensorExpr(a,b));
}

/* B(j)*A(j,i) -> Tensor1 */

template<class A, class B, class T, class U, int Dim0, int Dim1, char i,char j>
inline const Tensor1_Expr<const Tensor2_times_Tensor1_0<A,B,T,U,Dim0,Dim1,i,j>,
  typename promote<T,U>::V,Dim0,i>
operator*(const Tensor1_Expr<B,U,Dim1,j> &b,
	  const Tensor2_Expr<A,T,Dim1,Dim0,j,i> &a)
{
  typedef const Tensor2_times_Tensor1_0<A,B,T,U,Dim0,Dim1,i,j> TensorExpr;
  return Tensor1_Expr<TensorExpr,typename promote<T,U>::V,Dim0,i>(TensorExpr(a,b));
}

/* A(i,j)*B(k) -> Tensor3 */

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k>
class Tensor2_times_Tensor1
{
  const Tensor2_Expr<A,T,Dim0,Dim1,i,j> iterA;
  const Tensor1_Expr<B,U,Dim2,k> iterB;
public:
  Tensor2_times_Tensor1(const Tensor2_Expr<A,T,Dim0,Dim1,i,j> &a,
			const Tensor1_Expr<B,U,Dim2,k> &b):
    iterA(a), iterB(b) {}
  typename promote<T,U>::V operator()(const int N1, const int N2, const int N3) const
  {
    return iterA(N1,N2)*iterB(N3);
  }
};

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k>
inline const Tensor3_Expr<const Tensor2_times_Tensor1<A,B,T,U,Dim0,Dim1,Dim2,i,j,k>,
  typename promote<T,U>::V,Dim0,Dim1,Dim2,i,j,k>
operator*(const Tensor2_Expr<A,T,Dim0,Dim1,i,j> &a,
	  const Tensor1_Expr<B,U,Dim2,k> &b)
{
  typedef const Tensor2_times_Tensor1<A,B,T,U,Dim0,Dim1,Dim2,i,j,k> TensorExpr;
  return Tensor3_Expr<TensorExpr,typename promote<T,U>::V,Dim0,Dim1,Dim2,i,j,k>
    (TensorExpr(a,b));
}

/* B(k)*A(i,j) -> Tensor3 */

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k>
inline const Tensor3_Expr<const Tensor2_times_Tensor1<A,B,T,U,Dim0,Dim1,Dim2,i,j,k>,
  typename promote<T,U>::V,Dim0,Dim1,Dim2,i,j,k>
operator*(const Tensor1_Expr<B,U,Dim2,k> &b,
	  const Tensor2_Expr<A,T,Dim0,Dim1,i,j> &a)
{
  typedef const Tensor2_times_Tensor1<A,B,T,U,Dim0,Dim1,Dim2,i,j,k> TensorExpr;
  return Tensor3_Expr<TensorExpr,typename promote<T,U>::V,Dim0,Dim1,Dim2,i,j,k>
    (TensorExpr(a,b));
}
