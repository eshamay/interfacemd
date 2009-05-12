#include <iostream>
#include "one_over_1_minus_x_fast.h"

int main()
{
  double y[]={0,1,2};
  double a1[]={2,3,4};
  double a2[]={5,6,7};
  double a3[]={8,9,10};
  double a4[]={11,12,13};
  double a5[]={14,15,16};
  double a6[]={17,18,19};
  double a7[]={20,21,22};
  double a8[]={23,24,25};
  double a9[]={26,27,28};

  for(int ii=0;ii<1000000;ii++)
    {
      func9(y,a1,a2,a3,a4,a5,a6,a7,a8,a9);
    }
  std::cout << y[0] << " " << y[1] << " " << y[2] << std::endl;
}
