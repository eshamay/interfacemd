/* Multiply a Tensor1 and a Tensor2 together but don't contract, yielding a
   Tensor2. */

/* A(i,j) & B(i) -> Tensor2 */

template<class A, class B, class T, class U, int Dim0, int Dim1, char i,char j>
class Tensor2_and_Tensor1_0
{
  const Tensor2_Expr<A,T,Dim0,Dim1,i,j> iterA;
  const Tensor1_Expr<B,U,Dim0,i> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return iterA(N1,N2)*iterB(N1);
  }

  Tensor2_and_Tensor1_0(const Tensor2_Expr<A,T,Dim0,Dim1,i,j> &a,
			const Tensor1_Expr<B,U,Dim0,i> &b):
    iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim0, int Dim1, char i,char j>
inline const Tensor2_Expr<const Tensor2_and_Tensor1_0<A,B,T,U,Dim0,Dim1,i,j>,
  typename promote<T,U>::V,Dim0,Dim1,i,j>
operator&(const Tensor2_Expr<A,T,Dim0,Dim1,i,j> &a,
	  const Tensor1_Expr<B,U,Dim0,i> &b)
{
  typedef const Tensor2_and_Tensor1_0<A,B,T,U,Dim0,Dim1,i,j> TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim0,Dim1,i,j>
    (TensorExpr(a,b));
}

/* B(i) & A(i,j) -> Tensor2 */

template<class A, class B, class T, class U, int Dim0, int Dim1, char i,char j>
inline const Tensor2_Expr<const Tensor2_and_Tensor1_0
<A,B,T,U,Dim0,Dim1,i,j>,typename promote<T,U>::V,Dim0,Dim1,i,j>
operator&(const Tensor1_Expr<B,U,Dim0,i> &b,
	  const Tensor2_Expr<A,T,Dim0,Dim1,i,j> &a)
{
  typedef const Tensor2_and_Tensor1_0<A,B,T,U,Dim0,Dim1,i,j> TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim0,Dim1,i,j>
    (TensorExpr(a,b));
}

/* A(i,j) & B(j) -> Tensor2 */

template<class A, class B, class T, class U, int Dim0, int Dim1,
  char i, char j>
class Tensor2_and_Tensor1_1
{
  const Tensor2_Expr<A,T,Dim0,Dim1,i,j> iterA;
  const Tensor1_Expr<B,U,Dim1,j> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return iterA(N1,N2)*iterB(N2);
  }

  Tensor2_and_Tensor1_1(const Tensor2_Expr<A,T,Dim0,Dim1,i,j> &a,
			const Tensor1_Expr<B,U,Dim1,j> &b):
    iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim0, int Dim1,
  char i, char j>
inline const Tensor2_Expr<const Tensor2_and_Tensor1_1
<A,B,T,U,Dim0,Dim1,i,j>,typename promote<T,U>::V,Dim0,Dim1,i,j>
operator&(const Tensor2_Expr<A,T,Dim0,Dim1,i,j> &a,
	  const Tensor1_Expr<B,U,Dim1,j> &b)
{
  typedef const Tensor2_and_Tensor1_1<A,B,T,U,Dim0,Dim1,i,j> TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim0,Dim1,i,j>
    (TensorExpr(a,b));
}

/* B(j) & A(i,j) -> Tensor2 */

template<class A, class B, class T, class U, int Dim0, int Dim1,
  char i, char j>
inline const Tensor2_Expr<const Tensor2_and_Tensor1_1<A,B,T,U,Dim0,Dim1,i,j>,
  typename promote<T,U>::V,Dim0,Dim1,i,j>
operator&(const Tensor1_Expr<B,U,Dim1,j> &b,
	  const Tensor2_Expr<A,T,Dim0,Dim1,i,j> &a)
{
  typedef const Tensor2_and_Tensor1_1<A,B,T,U,Dim0,Dim1,i,j> TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim0,Dim1,i,j>
    (TensorExpr(a,b));
}
