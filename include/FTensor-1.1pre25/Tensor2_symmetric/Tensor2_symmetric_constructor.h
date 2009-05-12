/* A helper class that allows simple initialization of the Tensor2_symmetric,
   but only if it has the correct number of elements. */

template<class T, int Tensor_Dim>
class Tensor2_symmetric_constructor;

template<class T>
class Tensor2_symmetric_constructor<T,2>
{
public:
  Tensor2_symmetric_constructor(T data[(2*3)/2], T d00, T d01, T d11)
  {
    data[0]=d00;
    data[1]=d01;
    data[2]=d11;
  }
};

template<class T>
class Tensor2_symmetric_constructor<T,3>
{
public:
  Tensor2_symmetric_constructor(T data[(3*4)/2], T d00, T d01, T d02,
				T d11, T d12, T d22)
  {
    data[0]=d00;
    data[1]=d01;
    data[2]=d02;
    data[3]=d11;
    data[4]=d12;
    data[5]=d22;
  }
};

template<class T>
class Tensor2_symmetric_constructor<T,4>
{
public:
  Tensor2_symmetric_constructor(T data[(4*5)/2], T d00, T d01, T d02, T d03,
				T d11, T d12, T d13, T d22, T d23, T d33)
  {
    data[0]=d00;
    data[1]=d01;
    data[2]=d02;
    data[3]=d03;
    data[4]=d11;
    data[5]=d12;
    data[6]=d13;
    data[7]=d22;
    data[8]=d23;
    data[9]=d33;
  }
};
