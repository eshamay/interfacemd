/* Creates a Tensor2_symmetric expression by contracting two Tensor2's
   together. There are different versions, depending on where the
   contracting indices are located (i.e. whether it is A(i,j)^B(j,k)
   or A(i,j)^B(k,j)).  The classes are numbered to differentiate
   between these.  Thus, A(i,j)^B(j,k) has 10 appended to the name
   because I count from 0. */

/* A(i,j)^B(j,k) */

template<class A, class B, class T, class U, int Dim, int Dim1,
  char i, char j, char k>
class Tensor2_carat_Tensor2_10
{
  const Tensor2_Expr<A,T,Dim,Dim1,i,j> iterA;
  const Tensor2_Expr<B,U,Dim1,Dim,j,k> iterB;

  template<int Current_Dim>
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<Current_Dim> &N) const
  {
    return iterA(N1,Current_Dim-1)*iterB(Current_Dim-1,N2)
      + eval(N1,N2,Number<Current_Dim-1>());
  }
  typename promote<T,U>::V eval(const int N1, const int N2, const Number<1> &N) const
  {
    return iterA(N1,0)*iterB(0,N2);
  }
public:
  Tensor2_carat_Tensor2_10(const Tensor2_Expr<A,T,Dim,Dim1,i,j> &a,
			   const Tensor2_Expr<B,U,Dim1,Dim,j,k> &b):
    iterA(a), iterB(b) {}
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return eval(N1,N2,Number<Dim1>());
  }
};

template<class A, class B, class T, class U, int Dim, int Dim1,
  char i, char j, char k>
inline const Tensor2_symmetric_Expr
<const Tensor2_carat_Tensor2_10<A,B,T,U,Dim,Dim1,i,j,k>,
  typename promote<T,U>::V,Dim,i,k>
operator^(const Tensor2_Expr<A,T,Dim,Dim1,i,j> &a,
	  const Tensor2_Expr<B,U,Dim1,Dim,j,k> &b)
{
  typedef const Tensor2_carat_Tensor2_10<A,B,T,U,Dim,Dim1,i,j,k> TensorExpr;
  return Tensor2_symmetric_Expr<TensorExpr,typename promote<T,U>::V,Dim,i,k>
    (TensorExpr(a,b));
}

/* A(i,j)^B(k,j) */

template<class A, class B, class T, class U, int Dim, int Dim1,
  char i, char j, char k>
class Tensor2_carat_Tensor2_11
{
  const Tensor2_Expr<A,T,Dim,Dim1,i,j> iterA;
  const Tensor2_Expr<B,U,Dim,Dim1,k,j> iterB;

  template<int Current_Dim>
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<Current_Dim> &N) const
  {
    return iterA(N1,Current_Dim-1)*iterB(N2,Current_Dim-1)
      + eval(N1,N2,Number<Current_Dim-1>());
  }
  typename promote<T,U>::V eval(const int N1, const int N2, const Number<1> &N) const
  {
    return iterA(N1,0)*iterB(N2,0);
  }
public:
  Tensor2_carat_Tensor2_11(const Tensor2_Expr<A,T,Dim,Dim1,i,j> &a,
			   const Tensor2_Expr<B,U,Dim,Dim1,k,j> &b):
    iterA(a), iterB(b) {}
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return eval(N1,N2,Number<Dim1>());
  }
};

template<class A, class B, class T, class U, int Dim, int Dim1,
  char i, char j, char k>
inline const Tensor2_symmetric_Expr
<const Tensor2_carat_Tensor2_11<A,B,T,U,Dim,Dim1,i,j,k>,
  typename promote<T,U>::V,Dim,i,k>
operator^(const Tensor2_Expr<A,T,Dim,Dim1,i,j> &a,
	  const Tensor2_Expr<B,U,Dim,Dim1,k,j> &b)
{
  typedef const Tensor2_carat_Tensor2_11<A,B,T,U,Dim,Dim1,i,j,k> TensorExpr;
  return Tensor2_symmetric_Expr<TensorExpr,typename promote<T,U>::V,Dim,i,k>
    (TensorExpr(a,b));
}

/* A(j,i)^B(j,k) */

template<class A, class B, class T, class U, int Dim, int Dim1,
  char i, char j, char k>
class Tensor2_carat_Tensor2_00
{
  const Tensor2_Expr<A,T,Dim1,Dim,j,i> iterA;
  const Tensor2_Expr<B,U,Dim1,Dim,j,k> iterB;

  template<int Current_Dim>
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<Current_Dim> &N) const
  {
    return iterA(Current_Dim-1,N1)*iterB(Current_Dim-1,N2)
      + eval(N1,N2,Number<Current_Dim-1>());
  }
  typename promote<T,U>::V eval(const int N1, const int N2, const Number<1> &N) const
  {
    return iterA(0,N1)*iterB(0,N2);
  }
public:
  Tensor2_carat_Tensor2_00(const Tensor2_Expr<A,T,Dim1,Dim,j,i> &a,
			   const Tensor2_Expr<B,U,Dim1,Dim,j,k> &b):
    iterA(a), iterB(b) {}
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return eval(N1,N2,Number<Dim1>());
  }
};

template<class A, class B, class T, class U, int Dim, int Dim1,
  char i, char j, char k>
inline const Tensor2_symmetric_Expr
<const Tensor2_carat_Tensor2_00<A,B,T,U,Dim,Dim1,i,j,k>,
  typename promote<T,U>::V,Dim,i,k>
operator^(const Tensor2_Expr<A,T,Dim1,Dim,j,i> &a,
	  const Tensor2_Expr<B,U,Dim1,Dim,j,k> &b)
{
  typedef const Tensor2_carat_Tensor2_00<A,B,T,U,Dim,Dim1,i,j,k> TensorExpr;
  return Tensor2_symmetric_Expr<TensorExpr,typename promote<T,U>::V,Dim,i,k>
    (TensorExpr(a,b));
}

/* A(j,i)^B(k,j) */

template<class A, class B, class T, class U, int Dim, int Dim1,
  char i, char j, char k>
class Tensor2_carat_Tensor2_01
{
  const Tensor2_Expr<A,T,Dim1,Dim,j,i> iterA;
  const Tensor2_Expr<B,U,Dim,Dim1,k,j> iterB;

  template<int Current_Dim>
  typename promote<T,U>::V eval(const int N1, const int N2,
		       const Number<Current_Dim> &N) const
  {
    return iterA(Current_Dim-1,N1)*iterB(N2,Current_Dim-1)
      + eval(N1,N2,Number<Current_Dim-1>());
  }
  typename promote<T,U>::V eval(const int N1, const int N2, const Number<1> &N) const
  {
    return iterA(0,N1)*iterB(N2,0);
  }
public:
  Tensor2_carat_Tensor2_01(const Tensor2_Expr<A,T,Dim1,Dim,j,i> &a,
			   const Tensor2_Expr<B,U,Dim,Dim1,k,j> &b):
    iterA(a), iterB(b) {}
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return eval(N1,N2,Number<Dim1>());
  }
};

template<class A, class B, class T, class U, int Dim, int Dim1,
  char i, char j, char k>
inline const Tensor2_symmetric_Expr
<const Tensor2_carat_Tensor2_01<A,B,T,U,Dim,Dim1,i,j,k>,
  typename promote<T,U>::V,Dim,i,k>
operator^(const Tensor2_Expr<A,T,Dim1,Dim,j,i> &a,
	  const Tensor2_Expr<B,U,Dim,Dim1,k,j> &b)
{
  typedef const Tensor2_carat_Tensor2_01<A,B,T,U,Dim,Dim1,i,j,k> TensorExpr;
  return Tensor2_symmetric_Expr<TensorExpr,typename promote<T,U>::V,Dim,i,k>
    (TensorExpr(a,b));
}

/* I don't think that this product actually gives a Tensor4_ddg. */

//  /* A(i,j)^B(k,l) -> Tensor4_ddg(i,k,j,l) */

//  template<class A, class B, class T, class U, int Dim, int Dim1,
//    char i, char j, char k>
//  class Tensor2_carat_Tensor2_0213
//  {
//    const Tensor2_Expr<A,T,Dim01,Dim23,i,j> iterA;
//    const Tensor2_Expr<B,U,Dim01,Dim23,k,l> iterB;
//  public:
//    Tensor2_carat_Tensor2_0213(const Tensor2_Expr<A,T,Dim01,Dim23,i,j> &a,
//  			     const Tensor2_Expr<B,U,Dim01,Dim23,k,l> &b):
//      iterA(a), iterB(b) {}
//    typename promote<T,U>::V operator()(const int N1, const int N2, const int N3,
//  			     const int N4) const
//    {
//      return iterA(N1,N3)*iterB(N2,N4);
//    }
//  };

//  template<class A, class B, class T, class U, int Dim01, int Dim23,
//    char i, char j, char k, char l>
//  inline const Tensor4_ddg_Expr<const Tensor2_carat_Tensor2_0213<A,B,T,U,Dim01,Dim23,i,j,k,l>,typename promote<T,U>::V,Dim01,Dim23,i,k,j,l>
//  operator^(const Tensor2_Expr<A,T,Dim01,Dim23,i,j> &a,
//  	  const Tensor2_Expr<B,U,Dim01,Dim23,k,l> &b)
//  {
//    typedef const Tensor2_carat_Tensor2_0213<A,B,T,U,Dim01,Dim23,i,j,k,l>
//      TensorExpr;
//    return Tensor4_ddg_Expr<TensorExpr,typename promote<T,U>::V,Dim01,Dim23,i,k,j,l>
//      (TensorExpr(a,b));
//  }
