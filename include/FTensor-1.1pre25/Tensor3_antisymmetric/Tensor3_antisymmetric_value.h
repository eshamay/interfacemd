/* A general version, not for pointers. */

template <class T, int Tensor_Dim0, int Tensor_Dim12>
class Tensor3_antisymmetric
{
  T data[Tensor_Dim0][(Tensor_Dim12*(Tensor_Dim12-1))/2];
public:
  Tensor3_antisymmetric() {}

  /* Tensor_Dim0=2, Tensor_Dim12=2 */
  Tensor3_antisymmetric(T d001, T d101)
  {
    Tensor3_antisymmetric_constructor<T,2,2>(data,d001,d101);
  }

  /* Tensor_Dim0=3, Tensor_Dim12=3 */
  Tensor3_antisymmetric(T d001, T d002, T d012,	T d101, T d102, T d112,
			T d201, T d202, T d212)
  {
    Tensor3_antisymmetric_constructor<T,3,3>
      (data,d001,d002,d012,d101,d102,d112,d201,d202,d212);
  }

  /* Tensor_Dim0=4, Tensor_Dim12=4 */
  Tensor3_antisymmetric(T d001, T d002, T d003, T d012, T d013, T d023,
			T d101, T d102, T d103, T d112, T d113, T d123,
			T d201, T d202, T d203, T d212, T d213, T d223)
  {
    Tensor3_antisymmetric_constructor<T,4,4>
      (data,d001,d002,d003,d012,d013,d023,d101,d102,d103,d112,d113,d123,
       d201,d202,d203,d212,d213,d223);
  }

  /* There are two ways of accessing the values inside,
     unsafe(int,int,int) and operator(int,int,int).
     unsafe(int,int,int) will give you a wrong answer if you aren't
     careful.  The problem is that we only store the minimal set of
     components, but some have different signs.  We can't return the
     negative of a component, and assign something to it, because that
     would assign something to a temporary.  To get the correct answer
     if you don't want to change the value, just use
     operator(int,int,int). */

  T & unsafe(const int N1, const int N2, const int N3)
  {
#ifdef FTENSOR_DEBUG
    if(N1>=Tensor_Dim0 || N1<0
       || N2>=Tensor_Dim12 || N2<0 || N3>=Tensor_Dim12 || N3<0
       || N2>=N3)
      {
	std::cerr << "Bad index in Tensor3_antisymmetric<T,"
		  << Tensor_Dim0 << "," << Tensor_Dim12 << ">.operator("
		  << N1 << "," << N2 << "," << N3 << ")" << std::endl;
	abort();
      }
#endif
    return data[N1][N3-1+(N2*(2*(Tensor_Dim12-1)-N2-1))/2];
  }

  T operator()(const int N1, const int N2, const int N3) const
  {
#ifdef FTENSOR_DEBUG
    if(N1>=Tensor_Dim0 || N1<0
       || N2>=Tensor_Dim12 || N2<0 || N3>=Tensor_Dim12 || N3<0)
      {
	std::cerr << "Bad index in Tensor3_antisymmetric<T,"
		  << Tensor_Dim0 << "," << Tensor_Dim12 << ">.operator("
		  << N1 << "," << N2 << "," << N3 << ") const" << std::endl;
	abort();
      }
#endif
    return N2<N3 ? data[N1][N3-1+(N2*(2*(Tensor_Dim12-1)-N2-1))/2]
      : (N2>N3 ? -data[N1][N2-1+(N3*(2*(Tensor_Dim12-1)-N3-1))/2] : 0.0);
  }

  /* These operator()'s are the first part in constructing template
     expressions. */

  template<char i, char j, char k, int Dim0, int Dim12>
  Tensor3_antisymmetric_Expr<Tensor3_antisymmetric<T,Tensor_Dim0,Tensor_Dim12>,
    T,Dim0,Dim12,i,j,k>
  operator()(const Index<i,Dim0> index1, const Index<j,Dim12> index2,
	     const Index<k,Dim12> index3)
  {
    return Tensor3_antisymmetric_Expr<Tensor3_antisymmetric
      <T,Tensor_Dim0,Tensor_Dim12>,T,Dim0,Dim12,i,j,k>(*this);
  }

  template<char i, char j, char k, int Dim0, int Dim12>
  Tensor3_antisymmetric_Expr<const Tensor3_antisymmetric
  <T,Tensor_Dim0,Tensor_Dim12>,T,Dim0,Dim12,i,j,k>
  operator()(const Index<i,Dim0> index1, const Index<j,Dim12> index2,
	     const Index<k,Dim12> index3) const
  {
    return Tensor3_antisymmetric_Expr<const Tensor3_antisymmetric
      <T,Tensor_Dim0,Tensor_Dim12>,T,Dim0,Dim12,i,j,k>(*this);
  }
};
