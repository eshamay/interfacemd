#include <iostream>
#include "../FTensor.h"

using namespace FTensor;

int main()
{
  Tensor1<double,3> y(0,1,2);
  Tensor1<double,3> x(2,3,4);
  Tensor1<double,3> n(5,6,7);
  Tensor2 t2(1,2,3,4,5,6,7,8,9);

  const Index<'i',3> i;
  const Index<'j'> j;
  const Number<0> N0;
  const Number<1> N1;
  const Number<2> N2;

  for(int ii=0;ii<100000000;ii++)
    {
      y(i)+=(x(i)+n(i)+t2(i,N0))
//    	-(x(i)+n(i)+t2(i,N0))+(x(i)+n(i)+t2(i,N0))
//    	-(x(i)+n(i)+t2(i,N0))+(x(i)+n(i)+t2(i,N0))
//    	-(x(i)+n(i)+t2(i,N0))+(x(i)+n(i)+t2(i,N0))
//    	-(x(i)+n(i)+t2(i,N0))+(x(i)+n(i)+t2(i,N0))
//    	-(x(i)+n(i)+t2(i,N0))+(x(i)+n(i)+t2(i,N0))
//    	-(x(i)+n(i)+t2(i,N0))+(x(i)+n(i)+t2(i,N0))
//    	-(x(i)+n(i)+t2(i,N0))+(x(i)+n(i)+t2(i,N0))
//    	-(x(i)+n(i)+t2(i,N0))+(x(i)+n(i)+t2(i,N0))
//    	-(x(i)+n(i)+t2(i,N0))+(x(i)+n(i)+t2(i,N0))
//  -(x(i)+n(i)+t2(i,N0))+(x(i)+n(i)+t2(i,N0))
//  -(x(i)+n(i)+t2(i,N0))+(x(i)+n(i)+t2(i,N0))
//  -(x(i)+n(i)+t2(i,N0))+(x(i)+n(i)+t2(i,N0))
//  -(x(i)+n(i)+t2(i,N0))+(x(i)+n(i)+t2(i,N0))
	;
    }
  std::cout << y(0) << " " << y(1) << " " << y(2) << std::endl;
}

