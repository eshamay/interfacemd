/* Multiplies two Tensor1's together yielding a T (int, double, etc.)
   or a Tensor2. */

/* A(i)*B(i) -> T */

template<class A, class B, class T, class U, int Dim, char i, int Current_Dim>
inline const typename promote<T,U>::V T1_times_T1(const Tensor1_Expr<A,T,Dim,i> &a,
				      const Tensor1_Expr<B,U,Dim,i> &b,
				      const Number<Current_Dim> &N)
{
  return a(Current_Dim-1)*b(Current_Dim-1) +
    T1_times_T1(a,b,Number<Current_Dim-1>());
}

template<class A, class B, class T, class U, int Dim, char i>
inline const typename promote<T,U>::V T1_times_T1(const Tensor1_Expr<A,T,Dim,i> &a,
				      const Tensor1_Expr<B,U,Dim,i> &b,
				      const Number<1> &N)
{
  return a(0)*b(0);
}

template<class A, class B, class T, class U, int Dim, char i>
inline const typename promote<T,U>::V operator*(const Tensor1_Expr<A,T,Dim,i> &a,
				 const Tensor1_Expr<B,U,Dim,i> &b)
{
  return T1_times_T1(a,b,Number<Dim>());
}

/* A(i)*B(j) -> Tensor2 */

template<class A, class B, class T, class U, int Dim0, int Dim1, char i,char j>
class Tensor1_times_Tensor1
{
  const Tensor1_Expr<A,T,Dim0,i> iterA;
  const Tensor1_Expr<B,U,Dim1,j> iterB;
public:
  Tensor1_times_Tensor1(const Tensor1_Expr<A,T,Dim0,i> &a,
			const Tensor1_Expr<B,U,Dim1,j> &b):
    iterA(a), iterB(b) {}
  typename promote<T,U>::V operator()(const int N1, const int N2) const
  {
    return iterA(N1)*iterB(N2);
  }
};

template<class A, class B, class T, class U, int Dim0, int Dim1, char i,char j>
inline const Tensor2_Expr<const Tensor1_times_Tensor1<A,B,T,U,Dim0,Dim1,i,j>,
  typename promote<T,U>::V,Dim0,Dim1,i,j>
operator*(const Tensor1_Expr<A,T,Dim0,i> &a, const Tensor1_Expr<B,U,Dim1,j> &b)
{
  typedef const Tensor1_times_Tensor1<A,B,T,U,Dim0,Dim1,i,j> TensorExpr;
  return Tensor2_Expr<TensorExpr,typename promote<T,U>::V,Dim0,Dim1,i,j>
    (TensorExpr(a,b));
}
