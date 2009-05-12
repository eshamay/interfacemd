/* A helper class that allows simple initialization of the
   Tensor3_antisymmetric, but only if it has the correct number of
   elements. */

template<class T, int Tensor_Dim0, int Tensor_Dim12>
class Tensor3_antisymmetric_constructor;

template<class T>
class Tensor3_antisymmetric_constructor<T,2,2>
{
public:
  Tensor3_antisymmetric_constructor(T data[2][1], T d001, T d101)
  {
    data[0][0]=d001;
    data[1][0]=d101;
  }
};

template<class T>
class Tensor3_antisymmetric_constructor<T,3,3>
{
public:
  Tensor3_antisymmetric_constructor(T data[3][3], T d001, T d002, T d012,
				    T d101, T d102, T d112,
				    T d201, T d202, T d212)
  {
    data[0][0]=d001;
    data[1][0]=d101;
    data[2][0]=d201;

    data[0][1]=d002;
    data[1][1]=d102;
    data[2][1]=d202;

    data[0][2]=d012;
    data[1][2]=d112;
    data[2][2]=d212;
  }
};

template<class T>
class Tensor3_antisymmetric_constructor<T,4,4>
{
public:
  Tensor3_antisymmetric_constructor
  (T data[10][4], T d001, T d002, T d003, T d012, T d013, T d023,
   T d101, T d102, T d103, T d112, T d113, T d123,
   T d201, T d202, T d203, T d212, T d213, T d223)
  {
    data[0][0]=d001;
    data[0][1]=d002;
    data[0][2]=d003;
    data[0][3]=d012;
    data[0][4]=d013;
    data[0][5]=d023;

    data[1][0]=d001;
    data[1][1]=d002;
    data[1][2]=d003;
    data[1][3]=d012;
    data[1][4]=d013;
    data[1][5]=d023;

    data[2][0]=d001;
    data[2][1]=d002;
    data[2][2]=d003;
    data[2][3]=d012;
    data[2][4]=d013;
    data[2][5]=d023;

    data[3][0]=d001;
    data[3][1]=d002;
    data[3][2]=d003;
    data[3][3]=d012;
    data[3][4]=d013;
    data[3][5]=d023;
  }
};

