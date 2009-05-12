/* Adds a Tensor3 to a Tensor3, yielding a Tensor3_dg. */

/* A(i,j,k)+B(i,k,j)->Tensor3_dg */

template<class A, class B, class T, class U, int Dim0, int Dim12,
  char i, char j, char k>
class Tensor3_or_Tensor3_12
{
  const Tensor3_Expr<A,T,Dim0,Dim12,Dim12,i,j,k> iterA;
  const Tensor3_Expr<B,U,Dim0,Dim12,Dim12,i,k,j> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2, const int N3) const
  {
    return iterA(N3,N1,N2)+iterB(N3,N2,N1);
  }

  Tensor3_or_Tensor3_12(const Tensor3_Expr<A,T,Dim0,Dim12,Dim12,i,j,k> &a,
			const Tensor3_Expr<B,U,Dim0,Dim12,Dim12,i,k,j> &b):
    iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim0, int Dim12,
  char i, char j, char k>
inline const Tensor3_dg_Expr
<const Tensor3_or_Tensor3_12<A,B,T,U,Dim0,Dim12,i,j,k>,
  typename promote<T,U>::V,Dim0,Dim12,j,k,i>
operator||(const Tensor3_Expr<A,T,Dim0,Dim12,Dim12,i,j,k> &a,
	   const Tensor3_Expr<B,U,Dim0,Dim12,Dim12,i,k,j> &b)
{
  typedef const Tensor3_or_Tensor3_12<A,B,T,U,Dim0,Dim12,i,j,k> TensorExpr;
  return Tensor3_dg_Expr<TensorExpr,typename promote<T,U>::V,Dim0,Dim12,j,k,i>
    (TensorExpr(a,b));
}

/* A(i,j,k)+B(k,j,i)->Tensor3_dg */

template<class A, class B, class T, class U, int Dim02, int Dim1,
  char i, char j, char k>
class Tensor3_or_Tensor3_02
{
  const Tensor3_Expr<A,T,Dim02,Dim1,Dim02,i,j,k> iterA;
  const Tensor3_Expr<B,U,Dim02,Dim1,Dim02,k,j,i> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2, const int N3) const
  {
    return iterA(N1,N3,N2)+iterB(N2,N3,N1);
  }

  Tensor3_or_Tensor3_02(const Tensor3_Expr<A,T,Dim02,Dim1,Dim02,i,j,k> &a,
			const Tensor3_Expr<B,U,Dim02,Dim1,Dim02,k,j,i> &b):
    iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim02, int Dim1,
  char i, char j, char k>
inline const Tensor3_dg_Expr
<const Tensor3_or_Tensor3_02<A,B,T,U,Dim02,Dim1,i,j,k>,
  typename promote<T,U>::V,Dim02,Dim1,i,k,j>
operator||(const Tensor3_Expr<A,T,Dim02,Dim1,Dim02,i,j,k> &a,
	   const Tensor3_Expr<B,U,Dim02,Dim1,Dim02,k,j,i> &b)
{
  typedef const Tensor3_or_Tensor3_02<A,B,T,U,Dim02,Dim1,i,j,k> TensorExpr;
  return Tensor3_dg_Expr<TensorExpr,typename promote<T,U>::V,Dim02,Dim1,i,k,j>
    (TensorExpr(a,b));
}

