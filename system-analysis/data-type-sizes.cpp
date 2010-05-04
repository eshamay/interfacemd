#include <iostream>
#include <limits>

using std::cout;
using std::endl;
using std::numeric_limits;

int main()
{
  cout<<"The values for data type short ranges from: "
    <<numeric_limits<short>::min()<<" to "<<numeric_limits<short>::max();
  cout<<endl;

  cout<<"The values for the data type int ranges from: "
    <<numeric_limits<int>::min()<<" to "<<numeric_limits<int>::max();
  cout<<endl;

  cout<<"The values for the data type long ranges from: "
    <<numeric_limits<long>::min()<<" to "<<numeric_limits<long>::max();
  cout<<endl;

  cout<<"The values for the data type float ranges from: "
    <<numeric_limits<float>::min()<<" to "<<numeric_limits<float>::max();
  cout<<endl;

  cout<<"The values for the data type double ranges from: "
    <<numeric_limits<double>::min()<<" to "<<numeric_limits<double>::max();
  cout<<endl;

  cout<<"The values for the data type long double ranges from: "
    <<numeric_limits<long double>::min()<<" to "<<numeric_limits<long double>::max();
  cout<<endl;

  return 0;
} 
