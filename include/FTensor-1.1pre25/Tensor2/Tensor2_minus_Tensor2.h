/* Subtracts a Tensor2 from a Tensor2, yielding a Tensor2. */

/* A(i,j)-B(i,j) */

template<class A, class B, class T, class U, int Dim0, int Dim1, char i,char j>
class Tensor2_minus_Tensor2_01
{
  const Tensor2_Expr<A,T,Dim0,Dim1,i,j> iterA;
  const Tensor2_Expr<B,U,Dim0,Dim1,i,j> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return iterA(N1,N2)-iterB(N1,N2);
  }

  Tensor2_minus_Tensor2_01(const Tensor2_Expr<A,T,Dim0,Dim1,i,j> &a,
			   const Tensor2_Expr<B,U,Dim0,Dim1,i,j> &b):
    iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim0, int Dim1, char i,char j>
inline const Tensor2_Expr<const Tensor2_minus_Tensor2_01<A,B,T,U,Dim0,Dim1,i,j>,
  typename promote<T,U>::V,Dim0,Dim1,i,j>
operator-(const Tensor2_Expr<A,T,Dim0,Dim1,i,j> &a,
	  const Tensor2_Expr<B,U,Dim0,Dim1,i,j> &b)
{
  typedef const Tensor2_minus_Tensor2_01<A,B,T,U,Dim0,Dim1,i,j> TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim0,Dim1,i,j>
    (TensorExpr(a,b));
}

/* A(i,j)-B(j,i) */

template<class A, class B, class T, class U, int Dim0, int Dim1, char i,char j>
class Tensor2_minus_Tensor2_10
{
  const Tensor2_Expr<A,T,Dim0,Dim1,i,j> iterA;
  const Tensor2_Expr<B,U,Dim0,Dim1,j,i> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return iterA(N1,N2)-iterB(N2,N1);
  }

  Tensor2_minus_Tensor2_10(const Tensor2_Expr<A,T,Dim0,Dim1,i,j> &a,
			  const Tensor2_Expr<B,U,Dim0,Dim1,j,i> &b):
    iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim0, int Dim1,
  char i, char j>
inline const Tensor2_Expr<const Tensor2_minus_Tensor2_10<A,B,T,U,Dim0,Dim1,i,j>,
  typename promote<T,U>::V,Dim0,Dim1,i,j>
operator-(const Tensor2_Expr<A,T,Dim0,Dim1,i,j> &a,
	  const Tensor2_Expr<B,U,Dim0,Dim1,j,i> &b)
{
  typedef const Tensor2_minus_Tensor2_10<A,B,T,U,Dim0,Dim1,i,j> TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim0,Dim1,i,j>
    (TensorExpr(a,b));
}
