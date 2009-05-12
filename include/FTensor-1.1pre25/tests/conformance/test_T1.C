#include <iostream>
#include "../../FTensor.h"
#include "test_for_zero.h"
using namespace FTensor;
using namespace std;

void test_T1(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	     Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	     Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	     const Tensor2<double,3,3> &t2_3,
	     Tensor2_symmetric<double,3> &t2s_1,
	     const Tensor2_symmetric<double,3> &t2s_2,
	     const Tensor2_symmetric<double,3> &t2s_3)
{
  Index<'i',3> i;
  Index<'j',3> j;
  Index<'k',3> k;
  Index<'l',3> l;
  Index<'m',3> m;
  Index<'n',3> n;

  Number<0> N0;
  Number<1> N1;
  Number<2> N2;

  /* Tensor1 test cases. */

  t1_1(i)=2*t1_2(i);
  test_for_zero(t1_1(0)-2*t1_2(0),"T*T1(0)");
  test_for_zero(t1_1(1)-2*t1_2(1),"T*T1(1)");
  test_for_zero(t1_1(2)-2*t1_2(2),"T*T1(2)");
  t1_1(i)+=3*t1_2(i);
  test_for_zero(t1_1(0)-5*t1_2(0),"T1+=T1(0)");
  test_for_zero(t1_1(1)-5*t1_2(1),"T1+=T1(1)");
  test_for_zero(t1_1(2)-5*t1_2(2),"T1+=T1(2)");
  t1_1(i)-=t1_2(i)/2;
  test_for_zero(t1_1(0)-4.5*t1_2(0),"T1-=T1(0)");
  test_for_zero(t1_1(1)-4.5*t1_2(1),"T1-=T1(1)");
  test_for_zero(t1_1(2)-4.5*t1_2(2),"T1-=T1(2)");
  t1_1(i)*=3;
  test_for_zero(t1_1(0)-13.5*t1_2(0),"T1(0)*=T");
  test_for_zero(t1_1(1)-13.5*t1_2(1),"T1(1)*=T");
  test_for_zero(t1_1(2)-13.5*t1_2(2),"T1(2)*=T");
  t1_1(i)/=4.5;
  test_for_zero(t1_1(0)-3*t1_2(0),"T1(0)/=T");
  test_for_zero(t1_1(1)-3*t1_2(1),"T1(1)/=T");
  test_for_zero(t1_1(2)-3*t1_2(2),"T1(2)/=T");
  t1_1(i)+=10;
  test_for_zero(t1_1(0)-3*t1_2(0)-10,"T1(0)+=T");
  test_for_zero(t1_1(1)-3*t1_2(1)-10,"T1(1)+=T");
  test_for_zero(t1_1(2)-3*t1_2(2)-10,"T1(2)+=T");
  t1_1(i)-=7;
  test_for_zero(t1_1(0)-3*t1_2(0)-3,"T1(0)-=T");
  test_for_zero(t1_1(1)-3*t1_2(1)-3,"T1(1)-=T");
  test_for_zero(t1_1(2)-3*t1_2(2)-3,"T1(2)-=T");
  t1_1(i)=t1_2(i)+t1_1(i);
  test_for_zero(t1_1(0)-4*t1_2(0)-3,"T1+T1(0)");
  test_for_zero(t1_1(1)-4*t1_2(1)-3,"T1+T1(1)");
  test_for_zero(t1_1(2)-4*t1_2(2)-3,"T1+T1(2)");
  t1_1(i)=-t1_2(i)-t1_1(i);
  test_for_zero(t1_1(0)+5*t1_2(0)+3,"-T1-T1(0)");
  test_for_zero(t1_1(1)+5*t1_2(1)+3,"-T1-T1(1)");
  test_for_zero(t1_1(2)+5*t1_2(2)+3,"-T1-T1(2)");
  t1_1(i)=t1_2(i)+10;
  test_for_zero(t1_1(0)-t1_2(0)-10,"T1(0)+T");
  test_for_zero(t1_1(1)-t1_2(1)-10,"T1(1)+T");
  test_for_zero(t1_1(2)-t1_2(2)-10,"T1(2)+T");
  t1_1(i)=t1_2(i)-10;
  test_for_zero(t1_1(0)-t1_2(0)+10,"T1(0)-T");
  test_for_zero(t1_1(1)-t1_2(1)+10,"T1(1)-T");
  test_for_zero(t1_1(2)-t1_2(2)+10,"T1(2)-T");
  test_for_zero(t1_1(i)*t1_2(i)-(t1_2(0)*(t1_2(0)-10) + t1_2(1)*(t1_2(1)-10)
				 + t1_2(2)*(t1_2(2)-10)),"T1(i)*T1(i)");
  t1_1(i)=10-t1_2(i);
  test_for_zero(t1_1(0)+t1_2(0)-10,"T-T1(0)");
  test_for_zero(t1_1(1)+t1_2(1)-10,"T-T1(1)");
  test_for_zero(t1_1(2)+t1_2(2)-10,"T-T1(2)");
  t1_1(i)=t1_2(i);
  test_for_zero(t1_1(0)-t1_2(0),"T1=T1(0)");
  test_for_zero(t1_1(1)-t1_2(1),"T1=T1(1)");
  test_for_zero(t1_1(2)-t1_2(2),"T1=T1(2)");
  t1_1(i)=10;
  test_for_zero(t1_1(0)-10,"T1(0)=T");
  test_for_zero(t1_1(1)-10,"T1(1)=T");
  test_for_zero(t1_1(2)-10,"T1(2)=T");

  t1_1(i)=(t1_1(i)&t1_2(i));
  test_for_zero(t1_1(0)-10*t1_2(0),"T1&T1(0)");
  test_for_zero(t1_1(1)-10*t1_2(1),"T1&T1(1)");
  test_for_zero(t1_1(2)-10*t1_2(2),"T1&T1(2)");

  cout << endl;
}
