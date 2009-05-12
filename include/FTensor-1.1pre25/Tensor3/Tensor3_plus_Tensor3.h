/* Adds a Tensor3 to a Tensor3, yielding a Tensor3. */

/* A(i,j,k)+B(i,j,k)->Tensor3 */

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k>
class Tensor3_plus_Tensor3
{
  const Tensor3_Expr<A,T,Dim0,Dim1,Dim2,i,j,k> iterA;
  const Tensor3_Expr<B,U,Dim0,Dim1,Dim2,i,j,k> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2, const int N3) const
  {
    return iterA(N1,N2,N3)+iterB(N1,N2,N3);
  }

  Tensor3_plus_Tensor3(const Tensor3_Expr<A,T,Dim0,Dim1,Dim2,i,j,k> &a,
		       const Tensor3_Expr<B,U,Dim0,Dim1,Dim2,i,j,k> &b)
    : iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k>
inline const Tensor3_Expr
<const Tensor3_plus_Tensor3<A,B,T,U,Dim0,Dim1,Dim2,i,j,k>,
  typename promote<T,U>::V,Dim0,Dim1,Dim2,i,j,k>
operator+(const Tensor3_Expr<A,T,Dim0,Dim1,Dim2,i,j,k> &a,
	  const Tensor3_Expr<B,U,Dim0,Dim1,Dim2,i,j,k> &b)
{
  typedef const Tensor3_plus_Tensor3<A,B,T,U,Dim0,Dim1,Dim2,i,j,k> TensorExpr;
  return Tensor3_Expr<TensorExpr,typename promote<T,U>::V,Dim0,Dim1,Dim2,i,j,k>
    (TensorExpr(a,b));
}

/* A(i,j,k)+B(i,k,j)->Tensor3 */

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k>
class Tensor3_plus_Tensor3_21
{
  const Tensor3_Expr<A,T,Dim0,Dim1,Dim2,i,j,k> iterA;
  const Tensor3_Expr<B,U,Dim0,Dim2,Dim1,i,k,j> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2, const int N3) const
  {
    return iterA(N1,N2,N3)+iterB(N1,N3,N2);
  }

  Tensor3_plus_Tensor3_21(const Tensor3_Expr<A,T,Dim0,Dim1,Dim2,i,j,k> &a,
			  const Tensor3_Expr<B,U,Dim0,Dim2,Dim1,i,k,j> &b)
    :iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k>
inline const Tensor3_Expr
<const Tensor3_plus_Tensor3_21<A,B,T,U,Dim0,Dim1,Dim2,i,j,k>,
  typename promote<T,U>::V,Dim0,Dim1,Dim2,i,j,k>
operator+(const Tensor3_Expr<A,T,Dim0,Dim1,Dim2,i,j,k> &a,
	  const Tensor3_Expr<B,U,Dim0,Dim2,Dim1,i,k,j> &b)
{
  typedef const Tensor3_plus_Tensor3_21<A,B,T,U,Dim0,Dim1,Dim2,i,j,k>
    TensorExpr;
  return Tensor3_Expr<TensorExpr,typename promote<T,U>::V,Dim0,Dim1,Dim2,i,j,k>
    (TensorExpr(a,b));
}

/* A(i,j,k)+B(j,i,k)->Tensor3 */

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k>
class Tensor3_plus_Tensor3_10
{
  const Tensor3_Expr<A,T,Dim0,Dim1,Dim2,i,j,k> iterA;
  const Tensor3_Expr<B,U,Dim1,Dim0,Dim2,j,i,k> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2, const int N3) const
  {
    return iterA(N1,N2,N3)+iterB(N2,N1,N3);
  }

  Tensor3_plus_Tensor3_10(const Tensor3_Expr<A,T,Dim0,Dim1,Dim2,i,j,k> &a,
			  const Tensor3_Expr<B,U,Dim1,Dim0,Dim2,j,i,k> &b)
    :iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k>
inline const Tensor3_Expr
<const Tensor3_plus_Tensor3_10<A,B,T,U,Dim0,Dim1,Dim2,i,j,k>,
  typename promote<T,U>::V,Dim0,Dim1,Dim2,i,j,k>
operator+(const Tensor3_Expr<A,T,Dim0,Dim1,Dim2,i,j,k> &a,
	  const Tensor3_Expr<B,U,Dim1,Dim0,Dim2,j,i,k> &b)
{
  typedef const Tensor3_plus_Tensor3_10<A,B,T,U,Dim0,Dim1,Dim2,i,j,k>
    TensorExpr;
  return Tensor3_Expr<TensorExpr,typename promote<T,U>::V,Dim0,Dim1,Dim2,i,j,k>
    (TensorExpr(a,b));
}

/* A(i,j,k)+B(k,i,j)->Tensor3 */

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k>
class Tensor3_plus_Tensor3_120
{
  const Tensor3_Expr<A,T,Dim0,Dim1,Dim2,i,j,k> iterA;
  const Tensor3_Expr<B,U,Dim2,Dim0,Dim1,k,i,j> iterB;
public:
  typename promote<T,U>::V operator()(const int N1, const int N2, const int N3) const
  {
    return iterA(N1,N2,N3)+iterB(N3,N1,N2);
  }

  Tensor3_plus_Tensor3_120(const Tensor3_Expr<A,T,Dim0,Dim1,Dim2,i,j,k> &a,
			   const Tensor3_Expr<B,U,Dim2,Dim0,Dim1,k,i,j> &b)
    : iterA(a), iterB(b) {}
};

template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
  char i, char j, char k>
inline const Tensor3_Expr
<const Tensor3_plus_Tensor3_120<A,B,T,U,Dim0,Dim1,Dim2,i,j,k>,
  typename promote<T,U>::V,Dim0,Dim1,Dim2,i,j,k>
operator+(const Tensor3_Expr<A,T,Dim0,Dim1,Dim2,i,j,k> &a,
	  const Tensor3_Expr<B,U,Dim2,Dim0,Dim1,k,i,j> &b)
{
  typedef const Tensor3_plus_Tensor3_120<A,B,T,U,Dim0,Dim1,Dim2,i,j,k>
    TensorExpr;
  return Tensor3_Expr<TensorExpr,typename promote<T,U>::V,Dim0,Dim1,Dim2,i,j,k>
    (TensorExpr(a,b));
}
