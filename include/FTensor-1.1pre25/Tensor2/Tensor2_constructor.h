/* A helper class that allows simple initialization of the Tensor2,
   but only if it has the correct number of elements.  There are
   specializations for column-major and row-major orderings. */

template<class T, int Tensor_Dim0, int Tensor_Dim1, Layout layout>
class Tensor2_constructor;

/* Column major versions */

template<class T>
class Tensor2_constructor<T,2,2,column_major>
{
public:
  Tensor2_constructor(T data[2][2],T d00,T d01,T d10,T d11)
  {
    data[0][0]=d00;
    data[0][1]=d01;
    data[1][0]=d10;
    data[1][1]=d11;
  }
};

template<class T>
class Tensor2_constructor<T,3,2,column_major>
{
public:
  Tensor2_constructor(T data[3][2],T d00, T d01, T d10, T d11,
		      T d20, T d21)
  {
    data[0][0]=d00;
    data[0][1]=d01;
    data[1][0]=d10;
    data[1][1]=d11;
    data[2][0]=d20;
    data[2][1]=d21;
  }
};

template<class T>
class Tensor2_constructor<T,3,3,column_major>
{
public:
  Tensor2_constructor(T data[3][3],T d00, T d01, T d02, T d10, T d11, T d12,
		      T d20, T d21, T d22)
  {
    data[0][0]=d00;
    data[0][1]=d01;
    data[0][2]=d02;
    data[1][0]=d10;
    data[1][1]=d11;
    data[1][2]=d12;
    data[2][0]=d20;
    data[2][1]=d21;
    data[2][2]=d22;
  }
};

template<class T>
class Tensor2_constructor<T,4,4,column_major>
{
public:
  Tensor2_constructor(T data[4][4], T d00, T d01, T d02, T d03, T d10, T d11,
		      T d12, T d13, T d20, T d21, T d22, T d23, T d30, T d31,
		      T d32, T d33)
  {
    data[0][0]=d00;
    data[0][1]=d01;
    data[0][2]=d02;
    data[0][3]=d03;
    data[1][0]=d10;
    data[1][1]=d11;
    data[1][2]=d12;
    data[1][3]=d13;
    data[2][0]=d20;
    data[2][1]=d21;
    data[2][2]=d22;
    data[2][3]=d23;
    data[3][0]=d30;
    data[3][1]=d31;
    data[3][2]=d32;
    data[3][3]=d33;
  }
};

/* Row major versions */

template<class T>
class Tensor2_constructor<T,2,2,row_major>
{
public:
  Tensor2_constructor(T data[2][2],T d00,T d01,T d10,T d11)
  {
    data[0][0]=d00;
    data[1][0]=d01;
    data[0][1]=d10;
    data[1][1]=d11;
  }
};

template<class T>
class Tensor2_constructor<T,3,2,row_major>
{
public:
  Tensor2_constructor(T data[2][3],T d00, T d01, T d10, T d11,
		      T d20, T d21)
  {
    data[0][0]=d00;
    data[1][0]=d01;
    data[0][1]=d10;
    data[1][1]=d11;
    data[0][2]=d20;
    data[1][2]=d21;
  }
};

template<class T>
class Tensor2_constructor<T,3,3,row_major>
{
public:
  Tensor2_constructor(T data[3][3],T d00, T d01, T d02, T d10, T d11, T d12,
		      T d20, T d21, T d22)
  {
    data[0][0]=d00;
    data[1][0]=d01;
    data[2][0]=d02;
    data[0][1]=d10;
    data[1][1]=d11;
    data[2][1]=d12;
    data[0][2]=d20;
    data[1][2]=d21;
    data[2][2]=d22;
  }
};

template<class T>
class Tensor2_constructor<T,4,4,row_major>
{
public:
  Tensor2_constructor(T data[4][4], T d00, T d01, T d02, T d03, T d10, T d11,
		      T d12, T d13, T d20, T d21, T d22, T d23, T d30, T d31,
		      T d32, T d33)
  {
    data[0][0]=d00;
    data[1][0]=d01;
    data[2][0]=d02;
    data[3][0]=d03;
    data[0][1]=d10;
    data[1][1]=d11;
    data[2][1]=d12;
    data[3][1]=d13;
    data[0][2]=d20;
    data[1][2]=d21;
    data[2][2]=d22;
    data[3][2]=d23;
    data[0][3]=d30;
    data[1][3]=d31;
    data[2][3]=d32;
    data[3][3]=d33;
  }
};
