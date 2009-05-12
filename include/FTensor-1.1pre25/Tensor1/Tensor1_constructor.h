/* A helper class that allows simple initialization of the Tensor1,
   but only if it has the correct number of elements. */

template<class T, int Tensor_Dim>
class Tensor1_constructor;

template<class T>
class Tensor1_constructor<T,2>
{
public:
  Tensor1_constructor(T data[], T d0, T d1)
  {
    data[0]=d0;
    data[1]=d1;
  }
};

template<class T>
class Tensor1_constructor<T,3>
{
public:
  Tensor1_constructor(T data[], T d0, T d1, T d2)
  {
    data[0]=d0;
    data[1]=d1;
    data[2]=d2;
  }
};

template<class T>
class Tensor1_constructor<T,4>
{
public:
  Tensor1_constructor(T data[], T d0, T d1, T d2, T d3)
  {
    data[0]=d0;
    data[1]=d1;
    data[2]=d2;
    data[3]=d3;
  }
};
