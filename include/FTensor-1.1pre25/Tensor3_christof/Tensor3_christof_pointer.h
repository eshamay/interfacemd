/* A version for pointers */

template <class T, int Tensor_Dim0, int Tensor_Dim12>
class Tensor3_christof<T*,Tensor_Dim0,Tensor_Dim12>
{
  mutable T * restrict  data[Tensor_Dim0][(Tensor_Dim12*(Tensor_Dim12+1))/2];
public:
  Tensor3_christof() {}

  /* Tensor_Dim0=Tensor_Dim12=2 */
  Tensor3_christof(T *d000, T *d100, T *d001, T *d101, T *d011, T *d111)
  {
    Tensor3_christof_constructor<T* restrict,Tensor_Dim0,Tensor_Dim12>
      (data,d000,d100,d001,d101,d011,d111);
  }

  /* Tensor_Dim0=Tensor_Dim12=3 */
  Tensor3_christof(T *d000, T *d100, T *d200, T *d001, T *d101, T *d201,
		   T *d002, T *d102, T *d202, T *d011, T *d111, T *d211,
		   T *d012, T *d112, T *d212, T *d022, T *d122, T *d222)
  {
    Tensor3_christof_constructor<T* restrict,Tensor_Dim0,Tensor_Dim12>
      (data,d000,d100,d200,d001,d101,d201,d002,d102,d202,d011,d111,d211,
       d012,d112,d212,d022,d122,d222);
  }

  /* Tensor_Dim0=Tensor_Dim12=4 */
  Tensor3_christof(T *d000, T *d100, T *d200, T *d300,
		   T *d001, T *d101, T *d201, T *d301,
		   T *d002, T *d102, T *d202, T *d302,
		   T *d003, T *d103, T *d203, T *d303,
		   T *d011, T *d111, T *d211, T *d311,
		   T *d012, T *d112, T *d212, T *d312,
		   T *d013, T *d113, T *d213, T *d313,
		   T *d022, T *d122, T *d222, T *d322,
		   T *d023, T *d123, T *d223, T *d323,
		   T *d033, T *d133, T *d233, T *d333)
  {
    Tensor3_christof_constructor<T* restrict,Tensor_Dim0,Tensor_Dim12>
      (data,d000,d100,d200,d300,d001,d101,d201,d301,
       d002,d102,d202,d302,d003,d103,d203,d303,
       d011,d111,d211,d311,d012,d112,d212,d312,
       d013,d113,d213,d313,d022,d122,d222,d322,
       d023,d123,d223,d323,d033,d133,d233,d333);
  }

  /* There are two operator(int,int,int)'s, one for non-consts that lets you
     change the value, and one for consts that doesn't. */

  T & operator()(const int N1, const int N2, const int N3)
  {
#ifdef FTENSOR_DEBUG
    if(N1>=Tensor_Dim0 || N1<0
       || N2>=Tensor_Dim12 || N2<0 || N3>=Tensor_Dim12 || N3<0)
      {
	std::cerr << "Bad index in Tensor3_christof<T*,"
		  << Tensor_Dim0 << "," << Tensor_Dim12 << ">.operator("
		  << N1 << "," << N2 << "," << N3 << ")" << std::endl;
	abort();
      }
#endif
    return N2>N3 ? *data[N1][N2+(N3*(2*Tensor_Dim12-N3-1))/2]
      : *data[N1][N3+(N2*(2*Tensor_Dim12-N2-1))/2];
  }

  T operator()(const int N1, const int N2, const int N3) const
  {
#ifdef FTENSOR_DEBUG
    if(N1>=Tensor_Dim0 || N1<0
       || N2>=Tensor_Dim12 || N2<0 || N3>=Tensor_Dim12 || N3<0)
      {
	std::cerr << "Bad index in Tensor3_christof<T*,"
		  << Tensor_Dim0 << "," << Tensor_Dim12 << ">.operator("
		  << N1 << "," << N2 << "," << N3 << ") const" << std::endl;
	abort();
      }
#endif
    return N2>N3 ? *data[N1][N2+(N3*(2*Tensor_Dim12-N3-1))/2]
      : *data[N1][N3+(N2*(2*Tensor_Dim12-N2-1))/2];
  }

  T* ptr(const int N1, const int N2, const int N3) const
  {
#ifdef FTENSOR_DEBUG
    if(N1>=Tensor_Dim0 || N1<0
       || N2>=Tensor_Dim12 || N2<0 || N3>=Tensor_Dim12 || N3<0)
      {
	std::cerr << "Bad index in Tensor3_christof<T*,"
		  << Tensor_Dim0 << "," << Tensor_Dim12 << ">.ptr("
		  << N1 << "," << N2 << "," << N3 << ")" << std::endl;
	abort();
      }
#endif
    return N2>N3 ? data[N1][N2+(N3*(2*Tensor_Dim12-N3-1))/2]
      : data[N1][N3+(N2*(2*Tensor_Dim12-N2-1))/2];
  }

  /* These operator()'s are the first part in constructing template
     expressions.  I mix up the indices here so that it behaves like a
     Tensor3_dg.  That way I don't have to have a separate wrapper
     class Tensor3_christof_Expr, which simplifies things. */

  template<char i, char j, char k, int Dim0, int Dim12>
  Tensor3_dg_Expr<Tensor3_christof<T*,Tensor_Dim0,Tensor_Dim12>,
    T,Dim12,Dim0,i,j,k> operator()
    (const Index<k,Dim0> index1, const Index<i,Dim12> index2,
     const Index<j,Dim12> index3)
  {
    return Tensor3_dg_Expr<Tensor3_christof<T*,Tensor_Dim0,Tensor_Dim12>,
      T,Dim12,Dim0,i,j,k>(*this);
  }

  template<char i, char j, char k, int Dim0, int Dim12>
  Tensor3_dg_Expr<const Tensor3_christof<T*,Tensor_Dim0,Tensor_Dim12>,
    T,Dim12,Dim0,i,j,k>
  operator()(const Index<k,Dim0> index1, const Index<i,Dim12> index2,
	     const Index<j,Dim12> index3) const
  {
    return Tensor3_dg_Expr<const Tensor3_christof<T*,Tensor_Dim0,Tensor_Dim12>,
      T,Dim12,Dim0,i,j,k>(*this);
  }

  /* These operators are for internal contractions. */

  /* const versions */

  template<char i, char j, int Dim0, int Dim12>
  inline Tensor1_Expr<const Tensor3_contracted_12
  <const Tensor3_christof<T*,Tensor_Dim0,Tensor_Dim12>,T,Dim12>,T,Dim0,i>
  operator()(const Index<i,Dim0> index1, const Index<j,Dim12> index2,
	     const Index<j,Dim12> index3) const
  {
    typedef const Tensor3_contracted_12
      <const Tensor3_christof<T*,Tensor_Dim0,Tensor_Dim12>,T,Dim12>
      TensorExpr;
    return Tensor1_Expr<TensorExpr,T,Dim0,i>(TensorExpr(*this));
  }

  template<char i, char j, int Dim02, int Dim1>
  inline Tensor1_Expr<const Tensor3_contracted_02
  <const Tensor3_christof<T*,Tensor_Dim0,Tensor_Dim12>,T,Dim02>,T,Dim1,i>
  operator()(const Index<j,Dim02> index1, const Index<i,Dim1> index2,
	     const Index<j,Dim02> index3) const
  {
    typedef const Tensor3_contracted_02
      <const Tensor3_christof<T*,Tensor_Dim0,Tensor_Dim12>,T,Dim02>
      TensorExpr;
    return Tensor1_Expr<TensorExpr,T,Dim1,i>(TensorExpr(*this));
  }

  template<char i, char j, int Dim01, int Dim2>
  inline Tensor1_Expr<const Tensor3_contracted_01
  <const Tensor3_christof<T*,Tensor_Dim0,Tensor_Dim12>,T,Dim01>,T,Dim2,i>
  operator()(const Index<j,Dim01> index1, const Index<j,Dim01> index2,
	     const Index<i,Dim2> index3) const
  {
    typedef const Tensor3_contracted_01
      <const Tensor3_christof<T*,Tensor_Dim0,Tensor_Dim12>,T,Dim01>
      TensorExpr;
    return Tensor1_Expr<TensorExpr,T,Dim2,i>(TensorExpr(*this));
  }

  /* non-const versions */

  template<char i, char j, int Dim0, int Dim12>
  inline Tensor1_Expr<const Tensor3_contracted_12
  <const Tensor3_christof<T*,Tensor_Dim0,Tensor_Dim12>,T,Dim12>,T,Dim0,i>
  operator()(const Index<i,Dim0> index1, const Index<j,Dim12> index2,
	     const Index<j,Dim12> index3)
  {
    typedef const Tensor3_contracted_12
      <const Tensor3_christof<T*,Tensor_Dim0,Tensor_Dim12>,T,Dim12>
      TensorExpr;
    return Tensor1_Expr<TensorExpr,T,Dim0,i>(TensorExpr(*this));
  }

  template<char i, char j, int Dim02, int Dim1>
  inline Tensor1_Expr<const Tensor3_contracted_02
  <const Tensor3_christof<T*,Tensor_Dim0,Tensor_Dim12>,T,Dim02>,T,Dim1,i>
  operator()(const Index<j,Dim02> index1, const Index<i,Dim1> index2,
	     const Index<j,Dim02> index3)
  {
    typedef const Tensor3_contracted_02
      <const Tensor3_christof<T*,Tensor_Dim0,Tensor_Dim12>,T,Dim02>
      TensorExpr;
    return Tensor1_Expr<TensorExpr,T,Dim1,i>(TensorExpr(*this));
  }

  template<char i, char j, int Dim01, int Dim2>
  inline Tensor1_Expr<const Tensor3_contracted_01
  <const Tensor3_christof<T*,Tensor_Dim0,Tensor_Dim12>,T,Dim01>,T,Dim2,i>
  operator()(const Index<j,Dim01> index1, const Index<j,Dim01> index2,
	     const Index<i,Dim2> index3)
  {
    typedef const Tensor3_contracted_01
      <const Tensor3_christof<T*,Tensor_Dim0,Tensor_Dim12>,T,Dim01>
      TensorExpr;
    return Tensor1_Expr<TensorExpr,T,Dim2,i>(TensorExpr(*this));
  }

  /* This is for expressions where a number is used for one slot, and
     an index for the others, yielding a Tensor2_symmetric_Expr. */

  template<char i, char j, int N, int Dim12>
  Tensor2_symmetric_Expr<const Tensor3_christof_number_0
  <const Tensor3_christof<T*,Tensor_Dim0,Tensor_Dim12>,T,N>,T,Dim12,i,j>
  operator()(const Number<N> n1, const Index<i,Dim12> index1,
	     const Index<j,Dim12> index2) const
  {
    typedef const Tensor3_christof_number_0
      <const Tensor3_christof<T*,Tensor_Dim0,Tensor_Dim12>,T,N> TensorExpr;
    return Tensor2_symmetric_Expr<TensorExpr,T,Dim12,i,j>(TensorExpr(*this));
  }

  /* An int in one spot, Index for the others, yielding a Tensor2.  I
     can use the same structure for both, since Tensor3_christof is
     symmetric on the last two indices. */

  template<char i, char j, int Dim0, int Dim2>
  Tensor2_Expr<const Tensor3_christof_numeral_1
  <const Tensor3_christof<T*,Tensor_Dim0,Tensor_Dim12>,T>,T,Dim0,Dim2,i,j>
  operator()(const Index<i,Dim0> index1, const int N,
	     const Index<j,Dim2> index2) const
  {
    typedef const Tensor3_christof_numeral_1
      <const Tensor3_christof<T*,Tensor_Dim0,Tensor_Dim12>,T> TensorExpr;
    return Tensor2_Expr<TensorExpr,T,Dim0,Dim2,i,j>(TensorExpr(*this,N));
  }

  template<char i, char j, int Dim0, int Dim2>
  Tensor2_Expr<const Tensor3_christof_numeral_1
  <const Tensor3_christof<T*,Tensor_Dim0,Tensor_Dim12>,T>,T,Dim0,Dim2,i,j>
  operator()(const Index<i,Dim0> index1, const Index<j,Dim2> index2,
	     const int N) const
  {
    typedef const Tensor3_christof_numeral_1
      <const Tensor3_christof<T*,Tensor_Dim0,Tensor_Dim12>,T> TensorExpr;
    return Tensor2_Expr<TensorExpr,T,Dim0,Dim2,i,j>(TensorExpr(*this,N));
  }

  /* The ++ operator increments the pointer, not the number that the
     pointer points to.  This allows iterating over a grid. */

  const Tensor3_christof<T*,Tensor_Dim0,Tensor_Dim12> & operator++() const
  {
    for(int i=0;i<Tensor_Dim0;++i)
      for(int j=0;j<(Tensor_Dim12*(Tensor_Dim12+1))/2;++j)
	++data[i][j];
    return *this;
  }
};
