#include <iostream>
#include "one_over_1_minus_x.h"

int main()
{
  Tensor1<double,3> y(0,1,2);
  Tensor1<double,3> a1(2,3,4);

  for(int ii=0;ii<100000000;ii++)
    {
          func1(y,a1);
    }
  std::cout << y(0) << " " << y(1) << " " << y(2) << std::endl;
}
