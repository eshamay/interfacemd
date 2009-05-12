/* Adds Tensor4+Tensor4 -> Tensor4 */

/* A(i,j,k,l)+B(i,j,k,l) */

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  int Dim3, char i, char j, char k, char l>
class Tensor4_plus_Tensor4
{
  const Tensor4_Expr<A,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> iterA;
  const Tensor4_Expr<B,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2, const int N3,
		    const int N4) const
  {
    return iterA(N1,N2,N3,N4)+iterB(N1,N2,N3,N4);
  }

  Tensor4_plus_Tensor4(const Tensor4_Expr<A,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &a,
		       const Tensor4_Expr<B,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &b)
    : iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  int Dim3, char i, char j, char k, char l>
inline Tensor4_Expr
<const Tensor4_plus_Tensor4<A,B,T,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l>,
  typename promote<T,U>::V,Dim0,Dim1,Dim2,Dim3,i,j,k,l>
operator+(const Tensor4_Expr<A,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &a,
	  const Tensor4_Expr<B,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &b)
{
  typedef const Tensor4_plus_Tensor4<A,B,T,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l>
    TensorExpr;
  return Tensor4_Expr<TensorExpr,typename promote<T,U>::V,Dim0,Dim1,Dim2,Dim3,i,j,k,l>
    (TensorExpr(a,b));
}

/* A(i,j,k,l)+B(l,k,i,j) */

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  int Dim3, char i, char j, char k, char l>
class Tensor4_plus_Tensor4_3201
{
  const Tensor4_Expr<A,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> iterA;
  const Tensor4_Expr<B,U,Dim2,Dim3,Dim0,Dim1,l,k,i,j> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2, const int N3,
		    const int N4) const
  {
    return iterA(N1,N2,N3,N4)+iterB(N4,N3,N1,N2);
  }

  Tensor4_plus_Tensor4_3201
  (const Tensor4_Expr<A,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &a,
   const Tensor4_Expr<B,U,Dim2,Dim3,Dim0,Dim1,l,k,i,j> &b)
    : iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  int Dim3, char i, char j, char k, char l>
inline Tensor4_Expr
<const Tensor4_plus_Tensor4_3201<A,B,T,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l>,
  typename promote<T,U>::V,Dim0,Dim1,Dim2,Dim3,i,j,k,l>
operator+(const Tensor4_Expr<A,T,Dim0,Dim1,Dim2,Dim3,i,j,k,l> &a,
	  const Tensor4_Expr<B,U,Dim2,Dim3,Dim0,Dim1,l,k,i,j> &b)
{
  typedef const Tensor4_plus_Tensor4_3201<A,B,T,U,Dim0,Dim1,Dim2,Dim3,i,j,k,l>
    TensorExpr;
  return Tensor4_Expr<TensorExpr,typename promote<T,U>::V,Dim0,Dim1,Dim2,Dim3,i,j,k,l>
    (TensorExpr(a,b));
}
