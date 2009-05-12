/* The general version, not for pointers. */

template <class T, int Tensor_Dim>
class Tensor1
{
  T data[Tensor_Dim];
public:
  /* Initializations for varying numbers of elements, with each one
     defined for a particular Tensor_Dim.  To initialize a different
     dimension, just add the appropriate constructor and call to
     the Tensor1_constructor constructor. */
  Tensor1(T d0, T d1)
  {
    Tensor1_constructor<T,Tensor_Dim>(data,d0,d1);
  }
  Tensor1(T d0, T d1, T d2)
  {
    Tensor1_constructor<T,Tensor_Dim>(data,d0,d1,d2);
  }
  Tensor1(T d0, T d1, T d2, T d3)
  {
    Tensor1_constructor<T,Tensor_Dim>(data,d0,d1,d2,d3);
  }

  Tensor1() {}

  /* There are two operator(int)'s, one for non-consts that lets you
     change the value, and one for consts that doesn't. */

  T & operator()(const int N)
  {
#ifdef FTENSOR_DEBUG
    if(N>=Tensor_Dim || N<0)
      {
	std::cerr << "Bad index in Tensor1<T," << Tensor_Dim
		  << ">.operator(" << N << ")" << std::endl;
	abort();
      }
#endif
    return data[N];
  }
  T operator()(const int N) const
  {
#ifdef FTENSOR_DEBUG
    if(N>=Tensor_Dim || N<0)
      {
	std::cerr << "Bad index in Tensor1<T," << Tensor_Dim
		  << ">.operator(" << N << ") const" << std::endl;
	abort();
      }
#endif
    return data[N];
  }

  /* These operator()'s are the first part in constructing template
     expressions.  They can be used to slice off lower dimensional
     parts. They are not entirely safe, since you can accidently use a
     higher dimension than what is really allowed (like Dim=5). */

  template<char i, int Dim>
  Tensor1_Expr<Tensor1<T,Tensor_Dim>,T,Dim,i>
  operator()(const Index<i,Dim> &index)
  {
    return Tensor1_Expr<Tensor1<T,Tensor_Dim>,T,Dim,i>(*this);
  }

  template<char i, int Dim>
  Tensor1_Expr<const Tensor1<T,Tensor_Dim>,T,Dim,i>
  operator()(const Index<i,Dim> &index) const
  {
    return Tensor1_Expr<const Tensor1<T,Tensor_Dim>,T,Dim,i>(*this);
  }
};
