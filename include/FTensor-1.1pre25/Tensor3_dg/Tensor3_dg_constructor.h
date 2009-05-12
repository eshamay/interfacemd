/* A helper class that allows simple initialization of the Tensor3_dg,
   but only if it has the correct number of elements. */

template<class T, int Tensor_Dim01, int Tensor_Dim2>
class Tensor3_dg_constructor;

template<class T>
class Tensor3_dg_constructor<T,2,2>
{
public:
  Tensor3_dg_constructor(T data[3][2], T d000, T d001, T d010, T d011,
			 T d110, T d111)
  {
    data[0][0]=d000;
    data[0][1]=d001;
    data[1][0]=d010;
    data[1][1]=d011;
    data[2][0]=d110;
    data[2][1]=d111;
  }
};

template<class T>
class Tensor3_dg_constructor<T,3,3>
{
public:
  Tensor3_dg_constructor(T data[6][3], T d000, T d001, T d002, T d010, T d011,
			 T d012, T d020, T d021, T d022, T d110, T d111,
			 T d112, T d120, T d121, T d122, T d220, T d221,
			 T d222)
  {
    data[0][0]=d000;
    data[1][0]=d010;
    data[2][0]=d020;
    data[3][0]=d110;
    data[4][0]=d120;
    data[5][0]=d220;

    data[0][1]=d001;
    data[1][1]=d011;
    data[2][1]=d021;
    data[3][1]=d111;
    data[4][1]=d121;
    data[5][1]=d221;

    data[0][2]=d002;
    data[1][2]=d012;
    data[2][2]=d022;
    data[3][2]=d112;
    data[4][2]=d122;
    data[5][2]=d222;
  }
};

template<class T>
class Tensor3_dg_constructor<T,4,4>
{
public:
  Tensor3_dg_constructor(T data[10][4], T d000, T d001, T d002, T d003,
			 T d010, T d011, T d012, T d013,
			 T d020, T d021, T d022, T d023,
			 T d030, T d031, T d032, T d033,
			 T d110, T d111, T d112, T d113,
			 T d120, T d121, T d122, T d123,
			 T d130, T d131, T d132, T d133,
			 T d220, T d221, T d222, T d223,
			 T d230, T d231, T d232, T d233,
			 T d330, T d331, T d332, T d333)
  {
    data[0][0]=d000;
    data[1][0]=d010;
    data[2][0]=d020;
    data[3][0]=d030;
    data[4][0]=d110;
    data[5][0]=d120;
    data[6][0]=d130;
    data[7][0]=d220;
    data[8][0]=d230;
    data[9][0]=d330;

    data[0][1]=d001;
    data[1][1]=d011;
    data[2][1]=d021;
    data[3][1]=d031;
    data[4][1]=d111;
    data[5][1]=d121;
    data[6][1]=d131;
    data[7][1]=d221;
    data[8][1]=d231;
    data[9][1]=d331;

    data[0][2]=d002;
    data[1][2]=d012;
    data[2][2]=d022;
    data[3][2]=d032;
    data[4][2]=d112;
    data[5][2]=d122;
    data[6][2]=d132;
    data[7][2]=d222;
    data[8][2]=d232;
    data[9][2]=d332;

    data[0][3]=d003;
    data[1][3]=d013;
    data[2][3]=d023;
    data[3][3]=d033;
    data[4][3]=d113;
    data[5][3]=d123;
    data[6][3]=d133;
    data[7][3]=d223;
    data[8][3]=d233;
    data[9][3]=d333;
  }
};

template<class T>
class Tensor3_dg_constructor<T,4,3>
{
public:
  Tensor3_dg_constructor(T data[10][3], T d000, T d001, T d002,
			 T d010, T d011, T d012,
			 T d020, T d021, T d022,
			 T d030, T d031, T d032,
			 T d110, T d111, T d112,
			 T d120, T d121, T d122,
			 T d130, T d131, T d132,
			 T d220, T d221, T d222,
			 T d230, T d231, T d232,
			 T d330, T d331, T d332)
  {
    data[0][0]=d000;
    data[1][0]=d010;
    data[2][0]=d020;
    data[3][0]=d030;
    data[4][0]=d110;
    data[5][0]=d120;
    data[6][0]=d130;
    data[7][0]=d220;
    data[8][0]=d230;
    data[9][0]=d330;

    data[0][1]=d001;
    data[1][1]=d011;
    data[2][1]=d021;
    data[3][1]=d031;
    data[4][1]=d111;
    data[5][1]=d121;
    data[6][1]=d131;
    data[7][1]=d221;
    data[8][1]=d231;
    data[9][1]=d331;

    data[0][2]=d002;
    data[1][2]=d012;
    data[2][2]=d022;
    data[3][2]=d032;
    data[4][2]=d112;
    data[5][2]=d122;
    data[6][2]=d132;
    data[7][2]=d222;
    data[8][2]=d232;
    data[9][2]=d332;
  }
};
