/* Subtracts a Tensor3_dg from a Tensor3_dg, yielding a Tensor3_dg or
   Tensor3. */

/* A(i,j,k)-B(i,j,k)->Tensor3_dg */

template<class A, class B, class T, class U, int Dim01, int Dim2,
  char i, char j, char k>
class Tensor3_dg_minus_Tensor3_dg
{
  const Tensor3_dg_Expr<A,T,Dim01,Dim2,i,j,k> iterA;
  const Tensor3_dg_Expr<B,U,Dim01,Dim2,i,j,k> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2, const int N3) const
  {
    return iterA(N1,N2,N3)-iterB(N1,N2,N3);
  }

  Tensor3_dg_minus_Tensor3_dg(const Tensor3_dg_Expr<A,T,Dim01,Dim2,i,j,k> &a, const Tensor3_dg_Expr<B,U,Dim01,Dim2,i,j,k> &b):
    iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim01, int Dim2,
  char i, char j, char k>
inline const Tensor3_dg_Expr<const Tensor3_dg_minus_Tensor3_dg
<A,B,T,U,Dim01,Dim2,i,j,k>,typename promote<T,U>::V,Dim01,Dim2,i,j,k>
operator-(const Tensor3_dg_Expr<A,T,Dim01,Dim2,i,j,k> &a,
	  const Tensor3_dg_Expr<B,U,Dim01,Dim2,i,j,k> &b)
{
  typedef const Tensor3_dg_minus_Tensor3_dg<A,B,T,U,Dim01,Dim2,i,j,k>
    TensorExpr;
  return Tensor3_dg_Expr<TensorExpr,typename promote<T,U>::V,Dim01,Dim2,i,j,k>
    (TensorExpr(a,b));
}

/* A(i,j,k)-B(k,j,i)->Tensor3 */

template<class A, class B, class T, class U, int Dim, char i, char j, char k>
class Tensor3_dg_minus_Tensor3_dg_02
{
  const Tensor3_dg_Expr<A,T,Dim,Dim,i,j,k> iterA;
  const Tensor3_dg_Expr<B,U,Dim,Dim,k,j,i> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2, const int N3) const
  {
    return iterA(N1,N2,N3)-iterB(N3,N2,N1);
  }

  Tensor3_dg_minus_Tensor3_dg_02(const Tensor3_dg_Expr<A,T,Dim,Dim,i,j,k> &a,
				const Tensor3_dg_Expr<B,U,Dim,Dim,k,j,i> &b):
    iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim, char i, char j, char k>
inline const Tensor3_Expr<const Tensor3_dg_minus_Tensor3_dg_02<A,B,T,U,Dim,i,j,k>,
  typename promote<T,U>::V,Dim,Dim,Dim,i,j,k>
operator-(const Tensor3_dg_Expr<A,T,Dim,Dim,i,j,k> &a,
	  const Tensor3_dg_Expr<B,U,Dim,Dim,k,j,i> &b)
{
  typedef const Tensor3_dg_minus_Tensor3_dg_02<A,B,T,U,Dim,i,j,k> TensorExpr;
  return Tensor3_Expr<TensorExpr,typename promote<T,U>::V,Dim,Dim,Dim,i,j,k>
    (TensorExpr(a,b));
}
