#include <iostream>
#include "one_over_1_minus_x.h"

int main()
{
  Tensor1<double,3> y(0,1,2);
  Tensor1<double,3> a1(2,3,4);
  Tensor1<double,3> a2(5,6,7);
  Tensor1<double,3> a3(8,9,10);
  Tensor1<double,3> a4(11,12,13);
  Tensor1<double,3> a5(14,15,16);
  Tensor1<double,3> a6(17,18,19);
  Tensor1<double,3> a7(20,21,22);
  Tensor1<double,3> a8(23,24,25);
  Tensor1<double,3> a9(26,27,28);

  for(int ii=0;ii<1000000;ii++)
    {
      func9(y,a1,a2,a3,a4,a5,a6,a7,a8,a9);
    }
  std::cout << y(0) << " " << y(1) << " " << y(2) << std::endl;
}
