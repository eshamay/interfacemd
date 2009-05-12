/* A version for pointers. */

template <class T, int Tensor_Dim01, int Tensor_Dim23>
class Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>
{
  mutable T * restrict data[(Tensor_Dim01*(Tensor_Dim01+1))/2][(Tensor_Dim23*(Tensor_Dim23+1))/2];
public:
  /* There are two operator(int,int,int,int)'s, one for non-consts
     that lets you change the value, and one for consts that
     doesn't. */

  T & operator()(const int N1, const int N2, const int N3, const int N4)
  {
#ifdef FTENSOR_DEBUG
    if(N1>=Tensor_Dim01 || N1<0 || N2>=Tensor_Dim01 || N2<0
       || N3>=Tensor_Dim23 || N3<0 || N4>=Tensor_Dim23 || N4<0)
      {
	std::cerr << "Bad index in Tensor3_dg<T*,"
		  << Tensor_Dim01 << "," << Tensor_Dim23 << ">.operator("
		  << N1 << "," << N2 << "," << N3 << "," << N4 << ")"
		  << std::endl;
	abort();
      }
#endif
    return N1>N2 ? (N3>N4 ? *data[N1+(N2*(2*Tensor_Dim01-N2-1))/2]
		    [N3+(N4*(2*Tensor_Dim23-N4-1))/2]
		    : *data[N1+(N2*(2*Tensor_Dim01-N2-1))/2]
		    [N4+(N3*(2*Tensor_Dim23-N3-1))/2])
      : (N3>N4 ? *data[N2+(N1*(2*Tensor_Dim01-N1-1))/2]
	 [N3+(N4*(2*Tensor_Dim23-N4-1))/2]
	 : *data[N2+(N1*(2*Tensor_Dim01-N1-1))/2]
	 [N4+(N3*(2*Tensor_Dim23-N3-1))/2]);
  }

  T operator()(const int N1, const int N2, const int N3, const int N4)
    const
  {
#ifdef FTENSOR_DEBUG
    if(N1>=Tensor_Dim01 || N1<0 || N2>=Tensor_Dim01 || N2<0
       || N3>=Tensor_Dim23 || N3<0 || N4>=Tensor_Dim23 || N4<0)
      {
	std::cerr << "Bad index in Tensor3_dg<T*,"
		  << Tensor_Dim01 << "," << Tensor_Dim23 << ">.operator("
		  << N1 << "," << N2 << "," << N3 << "," << N4 << ") const"
		  << std::endl;
	abort();
      }
#endif
    return N1>N2 ? (N3>N4 ? *data[N1+(N2*(2*Tensor_Dim01-N2-1))/2]
		    [N3+(N4*(2*Tensor_Dim23-N4-1))/2]
		    : *data[N1+(N2*(2*Tensor_Dim01-N2-1))/2]
		    [N4+(N3*(2*Tensor_Dim23-N3-1))/2])
      : (N3>N4 ? *data[N2+(N1*(2*Tensor_Dim01-N1-1))/2]
	 [N3+(N4*(2*Tensor_Dim23-N4-1))/2]
	 : *data[N2+(N1*(2*Tensor_Dim01-N1-1))/2]
	 [N4+(N3*(2*Tensor_Dim23-N3-1))/2]);
  }

  T* ptr(const int N1, const int N2, const int N3, const int N4) const
  {
#ifdef FTENSOR_DEBUG
    if(N1>=Tensor_Dim01 || N1<0 || N2>=Tensor_Dim01 || N2<0
       || N3>=Tensor_Dim23 || N3<0 || N4>=Tensor_Dim23 || N4<0)
      {
	std::cerr << "Bad index in Tensor3_dg<T,"
		  << Tensor_Dim01 << "," << Tensor_Dim23 << ">.ptr("
		  << N1 << "," << N2 << "," << N3 << "," << N4 << ")"
		  << std::endl;
	abort();
      }
#endif
    return N1>N2 ? (N3>N4 ? data[N1+(N2*(2*Tensor_Dim01-N2-1))/2]
		    [N3+(N4*(2*Tensor_Dim23-N4-1))/2]
		    : data[N1+(N2*(2*Tensor_Dim01-N2-1))/2]
		    [N4+(N3*(2*Tensor_Dim23-N3-1))/2])
      : (N3>N4 ? data[N2+(N1*(2*Tensor_Dim01-N1-1))/2]
	 [N3+(N4*(2*Tensor_Dim23-N4-1))/2]
	 : data[N2+(N1*(2*Tensor_Dim01-N1-1))/2]
	 [N4+(N3*(2*Tensor_Dim23-N3-1))/2]);
  }

  /* These operator()'s are the first part in constructing template
     expressions.  They can be used to slice off lower dimensional
     parts. They are not entirely safe, since you can accidently use a
     higher dimension than what is really allowed (like Dim=5). */

  template<char i, char j, char k, char l, int Dim01, int Dim23>
  Tensor4_ddg_Expr<Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>,
    T,Dim01,Dim23,i,j,k,l> operator()
    (const Index<i,Dim01> index1, const Index<j,Dim01> index2,
     const Index<k,Dim23> index3, const Index<l,Dim23> index4)
  {
    return Tensor4_ddg_Expr<Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>,
      T,Dim01,Dim23,i,j,k,l> (*this);
  }


  template<char i, char j, char k, char l, int Dim01, int Dim23>
  const Tensor4_ddg_Expr<const Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>,
    T,Dim01,Dim23,i,j,k,l> operator()
    (const Index<i,Dim01> index1, const Index<j,Dim01> index2,
     const Index<k,Dim23> index3, const Index<l,Dim23> index4) const
  {
    return Tensor4_ddg_Expr<const Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>,
      T,Dim01,Dim23,i,j,k,l> (*this);
  }

  /* This is for expressions where a number is used for two slots, and
     an index for the other two, yielding a Tensor2_symmetric_Expr. */

  template<char i, char j, int N0, int N1, int Dim>
  const Tensor2_symmetric_Expr<const Tensor4_ddg_number_01
  <const Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>,T,N0,N1>,T,Dim,i,j>
  operator()(const Number<N0> n1, const Number<N1> n2,
	     const Index<i,Dim> index1, const Index<j,Dim> index2) const
  {
    typedef const Tensor4_ddg_number_01<const Tensor4_ddg
      <T*,Tensor_Dim01,Tensor_Dim23>,T,N0,N1> TensorExpr;
    return Tensor2_symmetric_Expr<TensorExpr,T,Dim,i,j>(TensorExpr(*this));
  }

  template<char i, char j, int N0, int N1, int Dim>
  Tensor2_symmetric_Expr<Tensor4_ddg_number_rhs_01
  <Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>,T,N0,N1>,T,Dim,i,j>
  operator()(const Number<N0> n1, const Number<N1> n2,
	     const Index<i,Dim> index1, const Index<j,Dim> index2)
  {
    typedef Tensor4_ddg_number_rhs_01<Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>,
      T,N0,N1> TensorExpr;
    return Tensor2_symmetric_Expr<TensorExpr,T,Dim,i,j>(*this);
  }

  /* This is for expressions where a number is used for one slot, and
     an index for the other three, yielding a Tensor3_dg_Expr. */

  template<char i, char j, char k, int N0, int Dim1, int Dim23>
  const Tensor3_dg_Expr<const Tensor4_ddg_number_0
  <const Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>,T,N0>,T,Dim23,Dim1,i,j,k>
  operator()(const Number<N0> n1, const Index<k,Dim1> index3,
	     const Index<i,Dim23> index1, const Index<j,Dim23> index2) const
  {
    typedef const Tensor4_ddg_number_0<const Tensor4_ddg
      <T*,Tensor_Dim01,Tensor_Dim23>,T,N0> TensorExpr;
    return Tensor3_dg_Expr<TensorExpr,T,Dim23,Dim1,i,j,k>(TensorExpr(*this));
  }

  template<char i, char j, char k, int N0, int Dim1, int Dim23>
  Tensor3_dg_Expr<Tensor4_ddg_number_rhs_0
  <Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>,T,N0>,T,Dim23,Dim1,i,j,k>
  operator()(const Number<N0> n1, const Index<k,Dim1> index3,
	     const Index<i,Dim23> index1, const Index<j,Dim23> index2)
  {
    typedef Tensor4_ddg_number_rhs_0<Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>,
      T,N0> TensorExpr;
    return Tensor3_dg_Expr<TensorExpr,T,Dim23,Dim1,i,j,k>(*this);
  }

  /* This is for expressions where an int (not a Number) is used for
     two slots, and an index for the other two, yielding a
     Tensor2_symmetric_Expr. */

  template<char i, char j, int Dim>
  const Tensor2_symmetric_Expr<const Tensor4_ddg_numeral_01
  <const Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>,T>,T,Dim,i,j>
  operator()(const int N0, const int N1,
	     const Index<i,Dim> index1, const Index<j,Dim> index2) const
  {
    typedef const Tensor4_ddg_numeral_01<const Tensor4_ddg
      <T*,Tensor_Dim01,Tensor_Dim23>,T> TensorExpr;
    return Tensor2_symmetric_Expr<TensorExpr,T,Dim,i,j>
      (TensorExpr(*this,N0,N1));
  }

  template<char i, char j, int Dim>
  const Tensor2_symmetric_Expr<const Tensor4_ddg_numeral_23
  <const Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>,T>,T,Dim,i,j>
  operator()(const Index<i,Dim> index1, const Index<j,Dim> index2,
	     const int N2, const int N3) const
  {
    typedef const Tensor4_ddg_numeral_23<const Tensor4_ddg
      <T*,Tensor_Dim01,Tensor_Dim23>,T> TensorExpr;
    return Tensor2_symmetric_Expr<TensorExpr,T,Dim,i,j>
      (TensorExpr(*this,N2,N3));
  }

  /* int in three slots, Index in the other yielding a Tensor1_Expr. */

  template<char i, int Dim>
  const Tensor1_Expr<const Tensor4_ddg_numeral_123<const Tensor4_ddg
  <T*,Tensor_Dim01,Tensor_Dim23>,T>,T,Dim,i>
  operator()(const Index<i,Dim> index1, const int N1, const int N2,
	     const int N3)
  {
    typedef const Tensor4_ddg_numeral_123<const Tensor4_ddg
      <T*,Tensor_Dim01,Tensor_Dim23>,T> TensorExpr;
    return Tensor1_Expr<TensorExpr,T,Dim,i>(TensorExpr(*this,N1,N2,N3));
  }

  template<char i, int Dim>
  const Tensor1_Expr<const Tensor4_ddg_numeral_123<const Tensor4_ddg
  <T*,Tensor_Dim01,Tensor_Dim23>,T>,T,Dim,i>
  operator()(const int N1, const Index<i,Dim> index1, const int N2,
	     const int N3)
  {
    typedef const Tensor4_ddg_numeral_123<const Tensor4_ddg
      <T*,Tensor_Dim01,Tensor_Dim23>,T> TensorExpr;
    return Tensor1_Expr<TensorExpr,T,Dim,i>(TensorExpr(*this,N1,N2,N3));
  }

  /* This is for expressions where an int (not a Number) is used for
     one slot, and an index for the other three, yielding a
     Tensor3_dg_Expr. */

  template<char i, char j, char k, int Dim1, int Dim23>
  const Tensor3_dg_Expr<const Tensor4_ddg_numeral_0
  <const Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23>,T>,T,Dim23,Dim1,i,j,k>
  operator()(const int N0, const Index<k,Dim1> index3,
	     const Index<i,Dim23> index1, const Index<j,Dim23> index2) const
  {
    typedef const Tensor4_ddg_numeral_0<const Tensor4_ddg
      <T*,Tensor_Dim01,Tensor_Dim23>,T> TensorExpr;
    return Tensor3_dg_Expr<TensorExpr,T,Dim23,Dim1,i,j,k>
      (TensorExpr(*this,N0));
  }

  /* The ++ operator increments the pointer, not the number that the
     pointer points to.  This allows iterating over a grid. */

  const Tensor4_ddg<T*,Tensor_Dim01,Tensor_Dim23> & operator++() const
  {
    for(int i=0;i<(Tensor_Dim01*(Tensor_Dim01+1))/2;++i)
      for(int j=0;j<(Tensor_Dim01*(Tensor_Dim01+1))/2;++j)
	++data[i][j];
    return *this;
  }
};
