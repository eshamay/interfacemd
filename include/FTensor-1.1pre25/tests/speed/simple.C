#include <iostream>
using namespace std;

void temp(int x, int y=10)
{
  cout << x << " " << y << endl;
}


//  template<class T, int Dim>
//  class helper;

//  template<class T>
//  class helper<T,1>
//  {
//  public:
//    helper(T temp[1], T a)
//    {
//      temp[0]=a;
//    }
//  };

//  template<class T, int Dim>
//  class simp
//  {
//  public:
//    T temp[Dim];
//    simp(T a)
//    {
//      helper<T,Dim>(temp,a);
//    }
//  };

int main()
{
  temp(10);
  temp(10,20);

//    simp<int,1> ttx(12);
//    cout << ttx.temp[0] << endl;
}
