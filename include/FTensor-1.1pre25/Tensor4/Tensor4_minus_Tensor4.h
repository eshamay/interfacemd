//  /* Subtracts a Tensor4 from a Tensor4, yielding a Tensor4_Riemann
//     (definitely not a general purpose subtraction routine). */

//  template<class A, class B, char i, char j, char k, char l>
//  class Tensor4_minus_Tensor4
//  {
//    Tensor4_Expr<A,i,j,k,l> iterA;
//    Tensor4_Expr<B,i,l,k,j> iterB;
//  public:
//    double operator()(const int N1, const int N2, const int N3,
//  			  const int N4) const
//    {
//      return iterA(N2,N3,N1,N4)-iterB(N2,N4,N1,N3);
//    }
//    Tensor4_minus_Tensor4(const Tensor4_Expr<A,i,j,k,l> &a,
//  			const Tensor4_Expr<B,i,l,k,j> &b):
//      iterA(a), iterB(b) {}
//  };

//  template<class A, class B, char i, char j, char k, char l>
//  inline Tensor4_Riemann_Expr<const Tensor4_minus_Tensor4
//  <const Tensor4_Expr<A,i,j,k,l>, const Tensor4_Expr<B,i,l,k,j>,i,j,k,l>,k,i,j,l>
//  operator-(const Tensor4_Expr<A,i,j,k,l> &a, const Tensor4_Expr<B,i,l,k,j> &b)
//  {
//    typedef const Tensor4_minus_Tensor4<const Tensor4_Expr<A,i,j,k,l>,
//      const Tensor4_Expr<B,i,l,k,j>,i,j,k,l> TensorExpr;
//    return Tensor4_Riemann_Expr<TensorExpr,k,i,j,l>(TensorExpr(a,b));
//  }
