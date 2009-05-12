#include <iostream>
#include "../../../FTensor.h"
#include "../test_for_zero.h"
using namespace FTensor;
using namespace std;

void test_T2VI(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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

  /* Tensor2 tests */

  t2_1(N0,i)*=4;
  t2_1(N1,i)*=6;
  t2_1(N2,i)*=8;
  test_for_zero(t2_1(0,0)-40,"T2(N,i)*=T(0,0)");
  test_for_zero(t2_1(0,1)-40,"T2(N,i)*=T(0,1)");
  test_for_zero(t2_1(0,2)-40,"T2(N,i)*=T(0,2)");
  test_for_zero(t2_1(1,0)-66,"T2(N,i)*=T(1,0)");
  test_for_zero(t2_1(1,1)-66,"T2(N,i)*=T(1,1)");
  test_for_zero(t2_1(1,2)-66,"T2(N,i)*=T(1,2)");
  test_for_zero(t2_1(2,0)-96,"T2(N,i)*=T(2,0)");
  test_for_zero(t2_1(2,1)-96,"T2(N,i)*=T(2,1)");
  test_for_zero(t2_1(2,2)-96,"T2(N,i)*=T(2,2)");
}
