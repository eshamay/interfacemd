/* A helper class that allows simple initialization of the Tensor3_christof,
   but only if it has the correct number of elements. */

template<class T, int Tensor_Dim0, int Tensor_Dim12>
class Tensor3_christof_constructor;

template<class T>
class Tensor3_christof_constructor<T,2,2>
{
public:
  Tensor3_christof_constructor(T data[2][3], T d000, T d100, T d001, T d101,
			       T d011, T d111)
  {
    data[0][0]=d000;
    data[1][0]=d100;
    data[0][1]=d001;
    data[1][1]=d101;
    data[0][2]=d011;
    data[1][2]=d111;
  }
};

template<class T>
class Tensor3_christof_constructor<T,3,3>
{
public:
  Tensor3_christof_constructor(T data[3][6], T d000, T d100, T d200, T d001, T d101,
			 T d201, T d002, T d102, T d202, T d011, T d111,
			 T d211, T d012, T d112, T d212, T d022, T d122,
			 T d222)
  {
    data[0][0]=d000;
    data[0][1]=d001;
    data[0][2]=d002;
    data[0][3]=d011;
    data[0][4]=d012;
    data[0][5]=d022;

    data[1][0]=d100;
    data[1][1]=d101;
    data[1][2]=d102;
    data[1][3]=d111;
    data[1][4]=d112;
    data[1][5]=d122;

    data[2][0]=d200;
    data[2][1]=d201;
    data[2][2]=d202;
    data[2][3]=d211;
    data[2][4]=d212;
    data[2][5]=d222;
  }
};

template<class T>
class Tensor3_christof_constructor<T,4,4>
{
public:
  Tensor3_christof_constructor(T data[4][10], T d000, T d100, T d200, T d300,
			 T d001, T d101, T d201, T d301,
			 T d002, T d102, T d202, T d302,
			 T d003, T d103, T d203, T d303,
			 T d011, T d111, T d211, T d311,
			 T d012, T d112, T d212, T d312,
			 T d013, T d113, T d213, T d313,
			 T d022, T d122, T d222, T d322,
			 T d023, T d123, T d223, T d323,
			 T d033, T d133, T d233, T d333)
  {
    data[0][0]=d000;
    data[0][1]=d001;
    data[0][2]=d002;
    data[0][3]=d003;
    data[0][4]=d011;
    data[0][5]=d012;
    data[0][6]=d013;
    data[0][7]=d022;
    data[0][8]=d023;
    data[0][9]=d033;

    data[1][0]=d100;
    data[1][1]=d101;
    data[1][2]=d102;
    data[1][3]=d103;
    data[1][4]=d111;
    data[1][5]=d112;
    data[1][6]=d113;
    data[1][7]=d122;
    data[1][8]=d123;
    data[1][9]=d133;

    data[2][0]=d200;
    data[2][1]=d201;
    data[2][2]=d202;
    data[2][3]=d203;
    data[2][4]=d211;
    data[2][5]=d212;
    data[2][6]=d213;
    data[2][7]=d222;
    data[2][8]=d223;
    data[2][9]=d233;

    data[3][0]=d300;
    data[3][1]=d301;
    data[3][2]=d302;
    data[3][3]=d303;
    data[3][4]=d311;
    data[3][5]=d312;
    data[3][6]=d313;
    data[3][7]=d322;
    data[3][8]=d323;
    data[3][9]=d333;
  }
};

template<class T>
class Tensor3_christof_constructor<T,3,4>
{
public:
  Tensor3_christof_constructor(T data[3][10], T d000, T d100, T d200,
			 T d001, T d101, T d201,
			 T d002, T d102, T d202,
			 T d003, T d103, T d203,
			 T d011, T d111, T d211,
			 T d012, T d112, T d212,
			 T d013, T d113, T d213,
			 T d022, T d122, T d222,
			 T d023, T d123, T d223,
			 T d033, T d133, T d233)
  {
    data[0][0]=d000;
    data[0][1]=d001;
    data[0][2]=d002;
    data[0][3]=d003;
    data[0][4]=d011;
    data[0][5]=d012;
    data[0][6]=d013;
    data[0][7]=d022;
    data[0][8]=d023;
    data[0][9]=d033;

    data[1][0]=d100;
    data[1][1]=d101;
    data[1][2]=d102;
    data[1][3]=d103;
    data[1][4]=d111;
    data[1][5]=d112;
    data[1][6]=d113;
    data[1][7]=d122;
    data[1][8]=d123;
    data[1][9]=d133;

    data[2][0]=d200;
    data[2][1]=d201;
    data[2][2]=d202;
    data[2][3]=d203;
    data[2][4]=d211;
    data[2][5]=d212;
    data[2][6]=d213;
    data[2][7]=d222;
    data[2][8]=d223;
    data[2][9]=d233;
  }
};
