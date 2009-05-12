#include <iostream>
#include "one_over_1_minus_x_fast.h"

int main()
{
  double y[]={0,1,2};
  double a1[]={2,3,4};
  double a2[]={5,6,7};
  double a3[]={8,9,10};

  for(int ii=0;ii<10000000;ii++)
    {
      func3(y,a1,a2,a3);
    }
  std::cout << y[0] << " " << y[1] << " " << y[2] << std::endl;
}
