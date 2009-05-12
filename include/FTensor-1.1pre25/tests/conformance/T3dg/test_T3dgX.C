#include <iostream>
#include "../../../FTensor.h"
#include "../test_for_zero.h"
using namespace FTensor;
using namespace std;

void test_T3dgX(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	       Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	       Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	       const Tensor2<double,3,3> &t2_3,
	       Tensor2_symmetric<double,3> &t2s_1,
	       const Tensor2_symmetric<double,3> &t2s_2,
	       const Tensor2_symmetric<double,3> &t2s_3,
	       Tensor3_dg<double,3,3> &t3dg_1,
	       const Tensor3_dg<double,3,3> &t3dg_2,
	       const Tensor3_dg<double,3,3> &t3dg_3,
	       Tensor3_christof<double,3,3> &t3ch_1,
	       const Tensor3_christof<double,3,3> &t3ch_2,
	       const Tensor3_christof<double,3,3> &t3ch_3)
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

  /* Tensor3_dg tests */

  t3dg_1(i,j,k)=10;
  test_for_zero(t3dg_1(0,0,0) - 10
		,"T3dg=T(0,0,0)");
  test_for_zero(t3dg_1(0,0,1) - 10
		,"T3dg=T(0,0,1)");
  test_for_zero(t3dg_1(0,0,2) - 10
		,"T3dg=T(0,0,2)");
  test_for_zero(t3dg_1(0,1,0) - 10
		,"T3dg=T(0,1,0)");
  test_for_zero(t3dg_1(0,1,1) - 10
		,"T3dg=T(0,1,1)");
  test_for_zero(t3dg_1(0,1,2) - 10
		,"T3dg=T(0,1,2)");
  test_for_zero(t3dg_1(0,2,0) - 10
		,"T3dg=T(0,2,0)");
  test_for_zero(t3dg_1(0,2,1) - 10
		,"T3dg=T(0,2,1)");
  test_for_zero(t3dg_1(0,2,2) - 10
		,"T3dg=T(0,2,2)");
  test_for_zero(t3dg_1(1,0,0) - 10
		,"T3dg=T(1,0,0)");
  test_for_zero(t3dg_1(1,0,1) - 10
		,"T3dg=T(1,0,1)");
  test_for_zero(t3dg_1(1,0,2) - 10
		,"T3dg=T(1,0,2)");
  test_for_zero(t3dg_1(1,1,0) - 10
		,"T3dg=T(1,1,0)");
  test_for_zero(t3dg_1(1,1,1) - 10
		,"T3dg=T(1,1,1)");
  test_for_zero(t3dg_1(1,1,2) - 10
		,"T3dg=T(1,1,2)");
  test_for_zero(t3dg_1(1,2,0) - 10
		,"T3dg=T(1,2,0)");
  test_for_zero(t3dg_1(1,2,1) - 10
		,"T3dg=T(1,2,1)");
  test_for_zero(t3dg_1(1,2,2) - 10
		,"T3dg=T(1,2,2)");
  test_for_zero(t3dg_1(2,0,0) - 10
		,"T3dg=T(2,0,0)");
  test_for_zero(t3dg_1(2,0,1) - 10
		,"T3dg=T(2,0,1)");
  test_for_zero(t3dg_1(2,0,2) - 10
		,"T3dg=T(2,0,2)");
  test_for_zero(t3dg_1(2,1,0) - 10
		,"T3dg=T(2,1,0)");
  test_for_zero(t3dg_1(2,1,1) - 10
		,"T3dg=T(2,1,1)");
  test_for_zero(t3dg_1(2,1,2) - 10
		,"T3dg=T(2,1,2)");
  test_for_zero(t3dg_1(2,2,0) - 10
		,"T3dg=T(2,2,0)");
  test_for_zero(t3dg_1(2,2,1) - 10
		,"T3dg=T(2,2,1)");
  test_for_zero(t3dg_1(2,2,2) - 10
		,"T3dg=T(2,2,2)");
}

