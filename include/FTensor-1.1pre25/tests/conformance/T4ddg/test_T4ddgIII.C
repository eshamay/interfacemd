#include <iostream>
#include "../../../FTensor.h"
#include "../test_for_zero.h"
using namespace FTensor;
using namespace std;

void test_T4ddgIII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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

  Tensor4_ddg<double,3,3> t4ddg_1, t4ddg_2, t4ddg_3;

  t4ddg_2(i,k,j,l)=t2s_2(i,k)*t2s_3(j,l);
  t2_1(i,k)=t4ddg_2(i,j,k,l)*t2s_2(j,l);
  test_for_zero(t2_1(0,0) - (t4ddg_2(0,0,0,0)*t2s_2(0,0)
			     + t4ddg_2(0,0,0,1)*t2s_2(0,1)
			     + t4ddg_2(0,0,0,2)*t2s_2(0,2)
			     + t4ddg_2(0,1,0,0)*t2s_2(1,0)
			     + t4ddg_2(0,1,0,1)*t2s_2(1,1)
			     + t4ddg_2(0,1,0,2)*t2s_2(1,2)
			     + t4ddg_2(0,2,0,0)*t2s_2(2,0)
			     + t4ddg_2(0,2,0,1)*t2s_2(2,1)
			     + t4ddg_2(0,2,0,2)*t2s_2(2,2))
		,"T4ddg(i,j,k,l)*T2s(j,l)(0,0)");
  test_for_zero(t2_1(0,1) - (t4ddg_2(0,0,1,0)*t2s_2(0,0)
			     + t4ddg_2(0,0,1,1)*t2s_2(0,1)
			     + t4ddg_2(0,0,1,2)*t2s_2(0,2)
			     + t4ddg_2(0,1,1,0)*t2s_2(1,0)
			     + t4ddg_2(0,1,1,1)*t2s_2(1,1)
			     + t4ddg_2(0,1,1,2)*t2s_2(1,2)
			     + t4ddg_2(0,2,1,0)*t2s_2(2,0)
			     + t4ddg_2(0,2,1,1)*t2s_2(2,1)
			     + t4ddg_2(0,2,1,2)*t2s_2(2,2))
		,"T4ddg(i,j,k,l)*T2s(j,l)(0,1)");
  test_for_zero(t2_1(0,2) - (t4ddg_2(0,0,2,0)*t2s_2(0,0)
			     + t4ddg_2(0,0,2,1)*t2s_2(0,1)
			     + t4ddg_2(0,0,2,2)*t2s_2(0,2)
			     + t4ddg_2(0,1,2,0)*t2s_2(1,0)
			     + t4ddg_2(0,1,2,1)*t2s_2(1,1)
			     + t4ddg_2(0,1,2,2)*t2s_2(1,2)
			     + t4ddg_2(0,2,2,0)*t2s_2(2,0)
			     + t4ddg_2(0,2,2,1)*t2s_2(2,1)
			     + t4ddg_2(0,2,2,2)*t2s_2(2,2))
		,"T4ddg(i,j,k,l)*T2s(j,l)(0,2)");
  test_for_zero(t2_1(1,0) - (t4ddg_2(1,0,0,0)*t2s_2(0,0)
			     + t4ddg_2(1,0,0,1)*t2s_2(0,1)
			     + t4ddg_2(1,0,0,2)*t2s_2(0,2)
			     + t4ddg_2(1,1,0,0)*t2s_2(1,0)
			     + t4ddg_2(1,1,0,1)*t2s_2(1,1)
			     + t4ddg_2(1,1,0,2)*t2s_2(1,2)
			     + t4ddg_2(1,2,0,0)*t2s_2(2,0)
			     + t4ddg_2(1,2,0,1)*t2s_2(2,1)
			     + t4ddg_2(1,2,0,2)*t2s_2(2,2))
		,"T4ddg(i,j,k,l)*T2s(j,l)(1,0)");
  test_for_zero(t2_1(1,1) - (t4ddg_2(1,0,1,0)*t2s_2(0,0)
			     + t4ddg_2(1,0,1,1)*t2s_2(0,1)
			     + t4ddg_2(1,0,1,2)*t2s_2(0,2)
			     + t4ddg_2(1,1,1,0)*t2s_2(1,0)
			     + t4ddg_2(1,1,1,1)*t2s_2(1,1)
			     + t4ddg_2(1,1,1,2)*t2s_2(1,2)
			     + t4ddg_2(1,2,1,0)*t2s_2(2,0)
			     + t4ddg_2(1,2,1,1)*t2s_2(2,1)
			     + t4ddg_2(1,2,1,2)*t2s_2(2,2))
		,"T4ddg(i,j,k,l)*T2s(j,l)(1,1)");
  test_for_zero(t2_1(1,2) - (t4ddg_2(1,0,2,0)*t2s_2(0,0)
			     + t4ddg_2(1,0,2,1)*t2s_2(0,1)
			     + t4ddg_2(1,0,2,2)*t2s_2(0,2)
			     + t4ddg_2(1,1,2,0)*t2s_2(1,0)
			     + t4ddg_2(1,1,2,1)*t2s_2(1,1)
			     + t4ddg_2(1,1,2,2)*t2s_2(1,2)
			     + t4ddg_2(1,2,2,0)*t2s_2(2,0)
			     + t4ddg_2(1,2,2,1)*t2s_2(2,1)
			     + t4ddg_2(1,2,2,2)*t2s_2(2,2))
		,"T4ddg(i,j,k,l)*T2s(j,l)(1,2)");
  test_for_zero(t2_1(2,0) - (t4ddg_2(2,0,0,0)*t2s_2(0,0)
			     + t4ddg_2(2,0,0,1)*t2s_2(0,1)
			     + t4ddg_2(2,0,0,2)*t2s_2(0,2)
			     + t4ddg_2(2,1,0,0)*t2s_2(1,0)
			     + t4ddg_2(2,1,0,1)*t2s_2(1,1)
			     + t4ddg_2(2,1,0,2)*t2s_2(1,2)
			     + t4ddg_2(2,2,0,0)*t2s_2(2,0)
			     + t4ddg_2(2,2,0,1)*t2s_2(2,1)
			     + t4ddg_2(2,2,0,2)*t2s_2(2,2))
		,"T4ddg(i,j,k,l)*T2s(j,l)(2,0)");
  test_for_zero(t2_1(2,1) - (t4ddg_2(2,0,1,0)*t2s_2(0,0)
			     + t4ddg_2(2,0,1,1)*t2s_2(0,1)
			     + t4ddg_2(2,0,1,2)*t2s_2(0,2)
			     + t4ddg_2(2,1,1,0)*t2s_2(1,0)
			     + t4ddg_2(2,1,1,1)*t2s_2(1,1)
			     + t4ddg_2(2,1,1,2)*t2s_2(1,2)
			     + t4ddg_2(2,2,1,0)*t2s_2(2,0)
			     + t4ddg_2(2,2,1,1)*t2s_2(2,1)
			     + t4ddg_2(2,2,1,2)*t2s_2(2,2))
		,"T4ddg(i,j,k,l)*T2s(j,l)(2,1)");
  test_for_zero(t2_1(2,2) - (t4ddg_2(2,0,2,0)*t2s_2(0,0)
			     + t4ddg_2(2,0,2,1)*t2s_2(0,1)
			     + t4ddg_2(2,0,2,2)*t2s_2(0,2)
			     + t4ddg_2(2,1,2,0)*t2s_2(1,0)
			     + t4ddg_2(2,1,2,1)*t2s_2(1,1)
			     + t4ddg_2(2,1,2,2)*t2s_2(1,2)
			     + t4ddg_2(2,2,2,0)*t2s_2(2,0)
			     + t4ddg_2(2,2,2,1)*t2s_2(2,1)
			     + t4ddg_2(2,2,2,2)*t2s_2(2,2))
		,"T4ddg(i,j,k,l)*T2s(j,l)(2,2)");

  t2_1(i,k)=t2s_3(j,l)*t4ddg_2(i,j,k,l);
  test_for_zero(t2_1(0,0) - (t4ddg_2(0,0,0,0)*t2s_3(0,0)
			     + t4ddg_2(0,0,0,1)*t2s_3(0,1)
			     + t4ddg_2(0,0,0,2)*t2s_3(0,2)
			     + t4ddg_2(0,1,0,0)*t2s_3(1,0)
			     + t4ddg_2(0,1,0,1)*t2s_3(1,1)
			     + t4ddg_2(0,1,0,2)*t2s_3(1,2)
			     + t4ddg_2(0,2,0,0)*t2s_3(2,0)
			     + t4ddg_2(0,2,0,1)*t2s_3(2,1)
			     + t4ddg_2(0,2,0,2)*t2s_3(2,2))
		,"T2s(j,l)*T4ddg(i,j,k,l)(0,0)");
  test_for_zero(t2_1(0,1) - (t4ddg_2(0,0,1,0)*t2s_3(0,0)
			     + t4ddg_2(0,0,1,1)*t2s_3(0,1)
			     + t4ddg_2(0,0,1,2)*t2s_3(0,2)
			     + t4ddg_2(0,1,1,0)*t2s_3(1,0)
			     + t4ddg_2(0,1,1,1)*t2s_3(1,1)
			     + t4ddg_2(0,1,1,2)*t2s_3(1,2)
			     + t4ddg_2(0,2,1,0)*t2s_3(2,0)
			     + t4ddg_2(0,2,1,1)*t2s_3(2,1)
			     + t4ddg_2(0,2,1,2)*t2s_3(2,2))
		,"T2s(j,l)*T4ddg(i,j,k,l)(0,1)");
  test_for_zero(t2_1(0,2) - (t4ddg_2(0,0,2,0)*t2s_3(0,0)
			     + t4ddg_2(0,0,2,1)*t2s_3(0,1)
			     + t4ddg_2(0,0,2,2)*t2s_3(0,2)
			     + t4ddg_2(0,1,2,0)*t2s_3(1,0)
			     + t4ddg_2(0,1,2,1)*t2s_3(1,1)
			     + t4ddg_2(0,1,2,2)*t2s_3(1,2)
			     + t4ddg_2(0,2,2,0)*t2s_3(2,0)
			     + t4ddg_2(0,2,2,1)*t2s_3(2,1)
			     + t4ddg_2(0,2,2,2)*t2s_3(2,2))
		,"T2s(j,l)*T4ddg(i,j,k,l)(0,2)");
  test_for_zero(t2_1(1,0) - (t4ddg_2(1,0,0,0)*t2s_3(0,0)
			     + t4ddg_2(1,0,0,1)*t2s_3(0,1)
			     + t4ddg_2(1,0,0,2)*t2s_3(0,2)
			     + t4ddg_2(1,1,0,0)*t2s_3(1,0)
			     + t4ddg_2(1,1,0,1)*t2s_3(1,1)
			     + t4ddg_2(1,1,0,2)*t2s_3(1,2)
			     + t4ddg_2(1,2,0,0)*t2s_3(2,0)
			     + t4ddg_2(1,2,0,1)*t2s_3(2,1)
			     + t4ddg_2(1,2,0,2)*t2s_3(2,2))
		,"T2s(j,l)*T4ddg(i,j,k,l)(1,0)");
  test_for_zero(t2_1(1,1) - (t4ddg_2(1,0,1,0)*t2s_3(0,0)
			     + t4ddg_2(1,0,1,1)*t2s_3(0,1)
			     + t4ddg_2(1,0,1,2)*t2s_3(0,2)
			     + t4ddg_2(1,1,1,0)*t2s_3(1,0)
			     + t4ddg_2(1,1,1,1)*t2s_3(1,1)
			     + t4ddg_2(1,1,1,2)*t2s_3(1,2)
			     + t4ddg_2(1,2,1,0)*t2s_3(2,0)
			     + t4ddg_2(1,2,1,1)*t2s_3(2,1)
			     + t4ddg_2(1,2,1,2)*t2s_3(2,2))
		,"T2s(j,l)*T4ddg(i,j,k,l)(1,1)");
  test_for_zero(t2_1(1,2) - (t4ddg_2(1,0,2,0)*t2s_3(0,0)
			     + t4ddg_2(1,0,2,1)*t2s_3(0,1)
			     + t4ddg_2(1,0,2,2)*t2s_3(0,2)
			     + t4ddg_2(1,1,2,0)*t2s_3(1,0)
			     + t4ddg_2(1,1,2,1)*t2s_3(1,1)
			     + t4ddg_2(1,1,2,2)*t2s_3(1,2)
			     + t4ddg_2(1,2,2,0)*t2s_3(2,0)
			     + t4ddg_2(1,2,2,1)*t2s_3(2,1)
			     + t4ddg_2(1,2,2,2)*t2s_3(2,2))
		,"T2s(j,l)*T4ddg(i,j,k,l)(1,2)");
  test_for_zero(t2_1(2,0) - (t4ddg_2(2,0,0,0)*t2s_3(0,0)
			     + t4ddg_2(2,0,0,1)*t2s_3(0,1)
			     + t4ddg_2(2,0,0,2)*t2s_3(0,2)
			     + t4ddg_2(2,1,0,0)*t2s_3(1,0)
			     + t4ddg_2(2,1,0,1)*t2s_3(1,1)
			     + t4ddg_2(2,1,0,2)*t2s_3(1,2)
			     + t4ddg_2(2,2,0,0)*t2s_3(2,0)
			     + t4ddg_2(2,2,0,1)*t2s_3(2,1)
			     + t4ddg_2(2,2,0,2)*t2s_3(2,2))
		,"T2s(j,l)*T4ddg(i,j,k,l)(2,0)");
  test_for_zero(t2_1(2,1) - (t4ddg_2(2,0,1,0)*t2s_3(0,0)
			     + t4ddg_2(2,0,1,1)*t2s_3(0,1)
			     + t4ddg_2(2,0,1,2)*t2s_3(0,2)
			     + t4ddg_2(2,1,1,0)*t2s_3(1,0)
			     + t4ddg_2(2,1,1,1)*t2s_3(1,1)
			     + t4ddg_2(2,1,1,2)*t2s_3(1,2)
			     + t4ddg_2(2,2,1,0)*t2s_3(2,0)
			     + t4ddg_2(2,2,1,1)*t2s_3(2,1)
			     + t4ddg_2(2,2,1,2)*t2s_3(2,2))
		,"T2s(j,l)*T4ddg(i,j,k,l)(2,1)");
  test_for_zero(t2_1(2,2) - (t4ddg_2(2,0,2,0)*t2s_3(0,0)
			     + t4ddg_2(2,0,2,1)*t2s_3(0,1)
			     + t4ddg_2(2,0,2,2)*t2s_3(0,2)
			     + t4ddg_2(2,1,2,0)*t2s_3(1,0)
			     + t4ddg_2(2,1,2,1)*t2s_3(1,1)
			     + t4ddg_2(2,1,2,2)*t2s_3(1,2)
			     + t4ddg_2(2,2,2,0)*t2s_3(2,0)
			     + t4ddg_2(2,2,2,1)*t2s_3(2,1)
			     + t4ddg_2(2,2,2,2)*t2s_3(2,2))
		,"T2s(j,l)*T4ddg(i,j,k,l)(2,2)");

  t2s_1(j,l)=t4ddg_2(i,k,j,l)*t2s_2(i,k);
  test_for_zero(t2s_1(0,0) - (t4ddg_2(0,0,0,0)*t2s_2(0,0)
			      + t4ddg_2(0,1,0,0)*t2s_2(0,1)
			      + t4ddg_2(0,2,0,0)*t2s_2(0,2)
			      + t4ddg_2(1,0,0,0)*t2s_2(1,0)
			      + t4ddg_2(1,1,0,0)*t2s_2(1,1)
			      + t4ddg_2(1,2,0,0)*t2s_2(1,2)
			      + t4ddg_2(2,0,0,0)*t2s_2(2,0)
			      + t4ddg_2(2,1,0,0)*t2s_2(2,1)
			      + t4ddg_2(2,2,0,0)*t2s_2(2,2))
		,"T4ddg(i,j,k,l)*T2s(i,k)(0,0)");
  test_for_zero(t2s_1(0,1) - (t4ddg_2(0,0,0,1)*t2s_2(0,0)
			      + t4ddg_2(0,1,0,1)*t2s_2(0,1)
			      + t4ddg_2(0,2,0,1)*t2s_2(0,2)
			      + t4ddg_2(1,0,0,1)*t2s_2(1,0)
			      + t4ddg_2(1,1,0,1)*t2s_2(1,1)
			      + t4ddg_2(1,2,0,1)*t2s_2(1,2)
			      + t4ddg_2(2,0,0,1)*t2s_2(2,0)
			      + t4ddg_2(2,1,0,1)*t2s_2(2,1)
			      + t4ddg_2(2,2,0,1)*t2s_2(2,2))
		,"T4ddg(i,j,k,l)*T2s(i,k)(0,1)");
  test_for_zero(t2s_1(0,2) - (t4ddg_2(0,0,0,2)*t2s_2(0,0)
			      + t4ddg_2(0,1,0,2)*t2s_2(0,1)
			      + t4ddg_2(0,2,0,2)*t2s_2(0,2)
			      + t4ddg_2(1,0,0,2)*t2s_2(1,0)
			      + t4ddg_2(1,1,0,2)*t2s_2(1,1)
			      + t4ddg_2(1,2,0,2)*t2s_2(1,2)
			      + t4ddg_2(2,0,0,2)*t2s_2(2,0)
			      + t4ddg_2(2,1,0,2)*t2s_2(2,1)
			      + t4ddg_2(2,2,0,2)*t2s_2(2,2))
		,"T4ddg(i,j,k,l)*T2s(i,k)(0,2)");
  test_for_zero(t2s_1(1,0) - (t4ddg_2(0,0,1,0)*t2s_2(0,0)
			      + t4ddg_2(0,1,1,0)*t2s_2(0,1)
			      + t4ddg_2(0,2,1,0)*t2s_2(0,2)
			      + t4ddg_2(1,0,1,0)*t2s_2(1,0)
			      + t4ddg_2(1,1,1,0)*t2s_2(1,1)
			      + t4ddg_2(1,2,1,0)*t2s_2(1,2)
			      + t4ddg_2(2,0,1,0)*t2s_2(2,0)
			      + t4ddg_2(2,1,1,0)*t2s_2(2,1)
			      + t4ddg_2(2,2,1,0)*t2s_2(2,2))
		,"T4ddg(i,j,k,l)*T2s(i,k)(1,0)");
  test_for_zero(t2s_1(1,1) - (t4ddg_2(0,0,1,1)*t2s_2(0,0)
			      + t4ddg_2(0,1,1,1)*t2s_2(0,1)
			      + t4ddg_2(0,2,1,1)*t2s_2(0,2)
			      + t4ddg_2(1,0,1,1)*t2s_2(1,0)
			      + t4ddg_2(1,1,1,1)*t2s_2(1,1)
			      + t4ddg_2(1,2,1,1)*t2s_2(1,2)
			      + t4ddg_2(2,0,1,1)*t2s_2(2,0)
			      + t4ddg_2(2,1,1,1)*t2s_2(2,1)
			      + t4ddg_2(2,2,1,1)*t2s_2(2,2))
		,"T4ddg(i,j,k,l)*T2s(i,k)(1,1)");
  test_for_zero(t2s_1(1,2) - (t4ddg_2(0,0,1,2)*t2s_2(0,0)
			      + t4ddg_2(0,1,1,2)*t2s_2(0,1)
			      + t4ddg_2(0,2,1,2)*t2s_2(0,2)
			      + t4ddg_2(1,0,1,2)*t2s_2(1,0)
			      + t4ddg_2(1,1,1,2)*t2s_2(1,1)
			      + t4ddg_2(1,2,1,2)*t2s_2(1,2)
			      + t4ddg_2(2,0,1,2)*t2s_2(2,0)
			      + t4ddg_2(2,1,1,2)*t2s_2(2,1)
			      + t4ddg_2(2,2,1,2)*t2s_2(2,2))
		,"T4ddg(i,j,k,l)*T2s(i,k)(1,2)");
  test_for_zero(t2s_1(2,0) - (t4ddg_2(0,0,2,0)*t2s_2(0,0)
			      + t4ddg_2(0,1,2,0)*t2s_2(0,1)
			      + t4ddg_2(0,2,2,0)*t2s_2(0,2)
			      + t4ddg_2(1,0,2,0)*t2s_2(1,0)
			      + t4ddg_2(1,1,2,0)*t2s_2(1,1)
			      + t4ddg_2(1,2,2,0)*t2s_2(1,2)
			      + t4ddg_2(2,0,2,0)*t2s_2(2,0)
			      + t4ddg_2(2,1,2,0)*t2s_2(2,1)
			      + t4ddg_2(2,2,2,0)*t2s_2(2,2))
		,"T4ddg(i,j,k,l)*T2s(i,k)(2,0)");
  test_for_zero(t2s_1(2,1) - (t4ddg_2(0,0,2,1)*t2s_2(0,0)
			      + t4ddg_2(0,1,2,1)*t2s_2(0,1)
			      + t4ddg_2(0,2,2,1)*t2s_2(0,2)
			      + t4ddg_2(1,0,2,1)*t2s_2(1,0)
			      + t4ddg_2(1,1,2,1)*t2s_2(1,1)
			      + t4ddg_2(1,2,2,1)*t2s_2(1,2)
			      + t4ddg_2(2,0,2,1)*t2s_2(2,0)
			      + t4ddg_2(2,1,2,1)*t2s_2(2,1)
			      + t4ddg_2(2,2,2,1)*t2s_2(2,2))
		,"T4ddg(i,j,k,l)*T2s(i,k)(2,1)");
  test_for_zero(t2s_1(2,2) - (t4ddg_2(0,0,2,2)*t2s_2(0,0)
			      + t4ddg_2(0,1,2,2)*t2s_2(0,1)
			      + t4ddg_2(0,2,2,2)*t2s_2(0,2)
			      + t4ddg_2(1,0,2,2)*t2s_2(1,0)
			      + t4ddg_2(1,1,2,2)*t2s_2(1,1)
			      + t4ddg_2(1,2,2,2)*t2s_2(1,2)
			      + t4ddg_2(2,0,2,2)*t2s_2(2,0)
			      + t4ddg_2(2,1,2,2)*t2s_2(2,1)
			      + t4ddg_2(2,2,2,2)*t2s_2(2,2))
		,"T4ddg(i,j,k,l)*T2s(i,k)(2,2)");

  t2s_1(j,l)=t2s_3(i,k)*t4ddg_2(i,k,j,l);
  test_for_zero(t2s_1(0,0) - (t4ddg_2(0,0,0,0)*t2s_3(0,0)
			      + t4ddg_2(0,1,0,0)*t2s_3(0,1)
			      + t4ddg_2(0,2,0,0)*t2s_3(0,2)
			      + t4ddg_2(1,0,0,0)*t2s_3(1,0)
			      + t4ddg_2(1,1,0,0)*t2s_3(1,1)
			      + t4ddg_2(1,2,0,0)*t2s_3(1,2)
			      + t4ddg_2(2,0,0,0)*t2s_3(2,0)
			      + t4ddg_2(2,1,0,0)*t2s_3(2,1)
			      + t4ddg_2(2,2,0,0)*t2s_3(2,2))
		,"T2s(i,k)*T4ddg(i,j,k,l)(0,0)");
  test_for_zero(t2s_1(0,1) - (t4ddg_2(0,0,0,1)*t2s_3(0,0)
			      + t4ddg_2(0,1,0,1)*t2s_3(0,1)
			      + t4ddg_2(0,2,0,1)*t2s_3(0,2)
			      + t4ddg_2(1,0,0,1)*t2s_3(1,0)
			      + t4ddg_2(1,1,0,1)*t2s_3(1,1)
			      + t4ddg_2(1,2,0,1)*t2s_3(1,2)
			      + t4ddg_2(2,0,0,1)*t2s_3(2,0)
			      + t4ddg_2(2,1,0,1)*t2s_3(2,1)
			      + t4ddg_2(2,2,0,1)*t2s_3(2,2))
		,"T2s(i,k)*T4ddg(i,j,k,l)(0,1)");
  test_for_zero(t2s_1(0,2) - (t4ddg_2(0,0,0,2)*t2s_3(0,0)
			      + t4ddg_2(0,1,0,2)*t2s_3(0,1)
			      + t4ddg_2(0,2,0,2)*t2s_3(0,2)
			      + t4ddg_2(1,0,0,2)*t2s_3(1,0)
			      + t4ddg_2(1,1,0,2)*t2s_3(1,1)
			      + t4ddg_2(1,2,0,2)*t2s_3(1,2)
			      + t4ddg_2(2,0,0,2)*t2s_3(2,0)
			      + t4ddg_2(2,1,0,2)*t2s_3(2,1)
			      + t4ddg_2(2,2,0,2)*t2s_3(2,2))
		,"T2s(i,k)*T4ddg(i,j,k,l)(0,2)");
  test_for_zero(t2s_1(1,0) - (t4ddg_2(0,0,1,0)*t2s_3(0,0)
			      + t4ddg_2(0,1,1,0)*t2s_3(0,1)
			      + t4ddg_2(0,2,1,0)*t2s_3(0,2)
			      + t4ddg_2(1,0,1,0)*t2s_3(1,0)
			      + t4ddg_2(1,1,1,0)*t2s_3(1,1)
			      + t4ddg_2(1,2,1,0)*t2s_3(1,2)
			      + t4ddg_2(2,0,1,0)*t2s_3(2,0)
			      + t4ddg_2(2,1,1,0)*t2s_3(2,1)
			      + t4ddg_2(2,2,1,0)*t2s_3(2,2))
		,"T2s(i,k)*T4ddg(i,j,k,l)(1,0)");
  test_for_zero(t2s_1(1,1) - (t4ddg_2(0,0,1,1)*t2s_3(0,0)
			      + t4ddg_2(0,1,1,1)*t2s_3(0,1)
			      + t4ddg_2(0,2,1,1)*t2s_3(0,2)
			      + t4ddg_2(1,0,1,1)*t2s_3(1,0)
			      + t4ddg_2(1,1,1,1)*t2s_3(1,1)
			      + t4ddg_2(1,2,1,1)*t2s_3(1,2)
			      + t4ddg_2(2,0,1,1)*t2s_3(2,0)
			      + t4ddg_2(2,1,1,1)*t2s_3(2,1)
			      + t4ddg_2(2,2,1,1)*t2s_3(2,2))
		,"T2s(i,k)*T4ddg(i,j,k,l)(1,1)");
  test_for_zero(t2s_1(1,2) - (t4ddg_2(0,0,1,2)*t2s_3(0,0)
			      + t4ddg_2(0,1,1,2)*t2s_3(0,1)
			      + t4ddg_2(0,2,1,2)*t2s_3(0,2)
			      + t4ddg_2(1,0,1,2)*t2s_3(1,0)
			      + t4ddg_2(1,1,1,2)*t2s_3(1,1)
			      + t4ddg_2(1,2,1,2)*t2s_3(1,2)
			      + t4ddg_2(2,0,1,2)*t2s_3(2,0)
			      + t4ddg_2(2,1,1,2)*t2s_3(2,1)
			      + t4ddg_2(2,2,1,2)*t2s_3(2,2))
		,"T2s(i,k)*T4ddg(i,j,k,l)(1,2)");
  test_for_zero(t2s_1(2,0) - (t4ddg_2(0,0,2,0)*t2s_3(0,0)
			      + t4ddg_2(0,1,2,0)*t2s_3(0,1)
			      + t4ddg_2(0,2,2,0)*t2s_3(0,2)
			      + t4ddg_2(1,0,2,0)*t2s_3(1,0)
			      + t4ddg_2(1,1,2,0)*t2s_3(1,1)
			      + t4ddg_2(1,2,2,0)*t2s_3(1,2)
			      + t4ddg_2(2,0,2,0)*t2s_3(2,0)
			      + t4ddg_2(2,1,2,0)*t2s_3(2,1)
			      + t4ddg_2(2,2,2,0)*t2s_3(2,2))
		,"T2s(i,k)*T4ddg(i,j,k,l)(2,0)");
  test_for_zero(t2s_1(2,1) - (t4ddg_2(0,0,2,1)*t2s_3(0,0)
			      + t4ddg_2(0,1,2,1)*t2s_3(0,1)
			      + t4ddg_2(0,2,2,1)*t2s_3(0,2)
			      + t4ddg_2(1,0,2,1)*t2s_3(1,0)
			      + t4ddg_2(1,1,2,1)*t2s_3(1,1)
			      + t4ddg_2(1,2,2,1)*t2s_3(1,2)
			      + t4ddg_2(2,0,2,1)*t2s_3(2,0)
			      + t4ddg_2(2,1,2,1)*t2s_3(2,1)
			      + t4ddg_2(2,2,2,1)*t2s_3(2,2))
		,"T2s(i,k)*T4ddg(i,j,k,l)(2,1)");
  test_for_zero(t2s_1(2,2) - (t4ddg_2(0,0,2,2)*t2s_3(0,0)
			      + t4ddg_2(0,1,2,2)*t2s_3(0,1)
			      + t4ddg_2(0,2,2,2)*t2s_3(0,2)
			      + t4ddg_2(1,0,2,2)*t2s_3(1,0)
			      + t4ddg_2(1,1,2,2)*t2s_3(1,1)
			      + t4ddg_2(1,2,2,2)*t2s_3(1,2)
			      + t4ddg_2(2,0,2,2)*t2s_3(2,0)
			      + t4ddg_2(2,1,2,2)*t2s_3(2,1)
			      + t4ddg_2(2,2,2,2)*t2s_3(2,2))
		,"T2s(i,k)*T4ddg(i,j,k,l)(2,2)");

  t4ddg_1(i,j,k,l)=(t4ddg_2(i,j,k,l)&t2s_3(i,j));
  test_for_zero(t4ddg_1(0,0,0,0) - (t4ddg_2(0,0,0,0)*t2s_3(0,0))
		,"T4ddg(i,j,k,l)&T2s(i,j)(0,0,0,0)");
  test_for_zero(t4ddg_1(0,0,0,1) - (t4ddg_2(0,0,0,1)*t2s_3(0,0))
		,"T4ddg(i,j,k,l)&T2s(i,j)(0,0,0,1)");
  test_for_zero(t4ddg_1(0,0,0,2) - (t4ddg_2(0,0,0,2)*t2s_3(0,0))
		,"T4ddg(i,j,k,l)&T2s(i,j)(0,0,0,2)");
  test_for_zero(t4ddg_1(0,0,1,0) - (t4ddg_2(0,0,1,0)*t2s_3(0,0))
		,"T4ddg(i,j,k,l)&T2s(i,j)(0,0,1,0)");
  test_for_zero(t4ddg_1(0,0,1,1) - (t4ddg_2(0,0,1,1)*t2s_3(0,0))
		,"T4ddg(i,j,k,l)&T2s(i,j)(0,0,1,1)");
  test_for_zero(t4ddg_1(0,0,1,2) - (t4ddg_2(0,0,1,2)*t2s_3(0,0))
		,"T4ddg(i,j,k,l)&T2s(i,j)(0,0,1,2)");
  test_for_zero(t4ddg_1(0,0,2,0) - (t4ddg_2(0,0,2,0)*t2s_3(0,0))
		,"T4ddg(i,j,k,l)&T2s(i,j)(0,0,2,0)");
  test_for_zero(t4ddg_1(0,0,2,1) - (t4ddg_2(0,0,2,1)*t2s_3(0,0))
		,"T4ddg(i,j,k,l)&T2s(i,j)(0,0,2,1)");
  test_for_zero(t4ddg_1(0,0,2,2) - (t4ddg_2(0,0,2,2)*t2s_3(0,0))
		,"T4ddg(i,j,k,l)&T2s(i,j)(0,0,2,2)");
  test_for_zero(t4ddg_1(0,1,0,0) - (t4ddg_2(0,1,0,0)*t2s_3(0,1))
		,"T4ddg(i,j,k,l)&T2s(i,j)(0,1,0,0)");
  test_for_zero(t4ddg_1(0,1,0,1) - (t4ddg_2(0,1,0,1)*t2s_3(0,1))
		,"T4ddg(i,j,k,l)&T2s(i,j)(0,1,0,1)");
  test_for_zero(t4ddg_1(0,1,0,2) - (t4ddg_2(0,1,0,2)*t2s_3(0,1))
		,"T4ddg(i,j,k,l)&T2s(i,j)(0,1,0,2)");
  test_for_zero(t4ddg_1(0,1,1,0) - (t4ddg_2(0,1,1,0)*t2s_3(0,1))
		,"T4ddg(i,j,k,l)&T2s(i,j)(0,1,1,0)");
  test_for_zero(t4ddg_1(0,1,1,1) - (t4ddg_2(0,1,1,1)*t2s_3(0,1))
		,"T4ddg(i,j,k,l)&T2s(i,j)(0,1,1,1)");
  test_for_zero(t4ddg_1(0,1,1,2) - (t4ddg_2(0,1,1,2)*t2s_3(0,1))
		,"T4ddg(i,j,k,l)&T2s(i,j)(0,1,1,2)");
  test_for_zero(t4ddg_1(0,1,2,0) - (t4ddg_2(0,1,2,0)*t2s_3(0,1))
		,"T4ddg(i,j,k,l)&T2s(i,j)(0,1,2,0)");
  test_for_zero(t4ddg_1(0,1,2,1) - (t4ddg_2(0,1,2,1)*t2s_3(0,1))
		,"T4ddg(i,j,k,l)&T2s(i,j)(0,1,2,1)");
  test_for_zero(t4ddg_1(0,1,2,2) - (t4ddg_2(0,1,2,2)*t2s_3(0,1))
		,"T4ddg(i,j,k,l)&T2s(i,j)(0,1,2,2)");
  test_for_zero(t4ddg_1(0,2,0,0) - (t4ddg_2(0,2,0,0)*t2s_3(0,2))
		,"T4ddg(i,j,k,l)&T2s(i,j)(0,2,0,0)");
  test_for_zero(t4ddg_1(0,2,0,1) - (t4ddg_2(0,2,0,1)*t2s_3(0,2))
		,"T4ddg(i,j,k,l)&T2s(i,j)(0,2,0,1)");
  test_for_zero(t4ddg_1(0,2,0,2) - (t4ddg_2(0,2,0,2)*t2s_3(0,2))
		,"T4ddg(i,j,k,l)&T2s(i,j)(0,2,0,2)");
  test_for_zero(t4ddg_1(0,2,1,0) - (t4ddg_2(0,2,1,0)*t2s_3(0,2))
		,"T4ddg(i,j,k,l)&T2s(i,j)(0,2,1,0)");
  test_for_zero(t4ddg_1(0,2,1,1) - (t4ddg_2(0,2,1,1)*t2s_3(0,2))
		,"T4ddg(i,j,k,l)&T2s(i,j)(0,2,1,1)");
  test_for_zero(t4ddg_1(0,2,1,2) - (t4ddg_2(0,2,1,2)*t2s_3(0,2))
		,"T4ddg(i,j,k,l)&T2s(i,j)(0,2,1,2)");
  test_for_zero(t4ddg_1(0,2,2,0) - (t4ddg_2(0,2,2,0)*t2s_3(0,2))
		,"T4ddg(i,j,k,l)&T2s(i,j)(0,2,2,0)");
  test_for_zero(t4ddg_1(0,2,2,1) - (t4ddg_2(0,2,2,1)*t2s_3(0,2))
		,"T4ddg(i,j,k,l)&T2s(i,j)(0,2,2,1)");
  test_for_zero(t4ddg_1(0,2,2,2) - (t4ddg_2(0,2,2,2)*t2s_3(0,2))
		,"T4ddg(i,j,k,l)&T2s(i,j)(0,2,2,2)");
  test_for_zero(t4ddg_1(1,0,0,0) - (t4ddg_2(1,0,0,0)*t2s_3(1,0))
		,"T4ddg(i,j,k,l)&T2s(i,j)(1,0,0,0)");
  test_for_zero(t4ddg_1(1,0,0,1) - (t4ddg_2(1,0,0,1)*t2s_3(1,0))
		,"T4ddg(i,j,k,l)&T2s(i,j)(1,0,0,1)");
  test_for_zero(t4ddg_1(1,0,0,2) - (t4ddg_2(1,0,0,2)*t2s_3(1,0))
		,"T4ddg(i,j,k,l)&T2s(i,j)(1,0,0,2)");
  test_for_zero(t4ddg_1(1,0,1,0) - (t4ddg_2(1,0,1,0)*t2s_3(1,0))
		,"T4ddg(i,j,k,l)&T2s(i,j)(1,0,1,0)");
  test_for_zero(t4ddg_1(1,0,1,1) - (t4ddg_2(1,0,1,1)*t2s_3(1,0))
		,"T4ddg(i,j,k,l)&T2s(i,j)(1,0,1,1)");
  test_for_zero(t4ddg_1(1,0,1,2) - (t4ddg_2(1,0,1,2)*t2s_3(1,0))
		,"T4ddg(i,j,k,l)&T2s(i,j)(1,0,1,2)");
  test_for_zero(t4ddg_1(1,0,2,0) - (t4ddg_2(1,0,2,0)*t2s_3(1,0))
		,"T4ddg(i,j,k,l)&T2s(i,j)(1,0,2,0)");
  test_for_zero(t4ddg_1(1,0,2,1) - (t4ddg_2(1,0,2,1)*t2s_3(1,0))
		,"T4ddg(i,j,k,l)&T2s(i,j)(1,0,2,1)");
  test_for_zero(t4ddg_1(1,0,2,2) - (t4ddg_2(1,0,2,2)*t2s_3(1,0))
		,"T4ddg(i,j,k,l)&T2s(i,j)(1,0,2,2)");
  test_for_zero(t4ddg_1(1,1,0,0) - (t4ddg_2(1,1,0,0)*t2s_3(1,1))
		,"T4ddg(i,j,k,l)&T2s(i,j)(1,1,0,0)");
  test_for_zero(t4ddg_1(1,1,0,1) - (t4ddg_2(1,1,0,1)*t2s_3(1,1))
		,"T4ddg(i,j,k,l)&T2s(i,j)(1,1,0,1)");
  test_for_zero(t4ddg_1(1,1,0,2) - (t4ddg_2(1,1,0,2)*t2s_3(1,1))
		,"T4ddg(i,j,k,l)&T2s(i,j)(1,1,0,2)");
  test_for_zero(t4ddg_1(1,1,1,0) - (t4ddg_2(1,1,1,0)*t2s_3(1,1))
		,"T4ddg(i,j,k,l)&T2s(i,j)(1,1,1,0)");
  test_for_zero(t4ddg_1(1,1,1,1) - (t4ddg_2(1,1,1,1)*t2s_3(1,1))
		,"T4ddg(i,j,k,l)&T2s(i,j)(1,1,1,1)");
  test_for_zero(t4ddg_1(1,1,1,2) - (t4ddg_2(1,1,1,2)*t2s_3(1,1))
		,"T4ddg(i,j,k,l)&T2s(i,j)(1,1,1,2)");
  test_for_zero(t4ddg_1(1,1,2,0) - (t4ddg_2(1,1,2,0)*t2s_3(1,1))
		,"T4ddg(i,j,k,l)&T2s(i,j)(1,1,2,0)");
  test_for_zero(t4ddg_1(1,1,2,1) - (t4ddg_2(1,1,2,1)*t2s_3(1,1))
		,"T4ddg(i,j,k,l)&T2s(i,j)(1,1,2,1)");
  test_for_zero(t4ddg_1(1,1,2,2) - (t4ddg_2(1,1,2,2)*t2s_3(1,1))
		,"T4ddg(i,j,k,l)&T2s(i,j)(1,1,2,2)");
  test_for_zero(t4ddg_1(1,2,0,0) - (t4ddg_2(1,2,0,0)*t2s_3(1,2))
		,"T4ddg(i,j,k,l)&T2s(i,j)(1,2,0,0)");
  test_for_zero(t4ddg_1(1,2,0,1) - (t4ddg_2(1,2,0,1)*t2s_3(1,2))
		,"T4ddg(i,j,k,l)&T2s(i,j)(1,2,0,1)");
  test_for_zero(t4ddg_1(1,2,0,2) - (t4ddg_2(1,2,0,2)*t2s_3(1,2))
		,"T4ddg(i,j,k,l)&T2s(i,j)(1,2,0,2)");
  test_for_zero(t4ddg_1(1,2,1,0) - (t4ddg_2(1,2,1,0)*t2s_3(1,2))
		,"T4ddg(i,j,k,l)&T2s(i,j)(1,2,1,0)");
  test_for_zero(t4ddg_1(1,2,1,1) - (t4ddg_2(1,2,1,1)*t2s_3(1,2))
		,"T4ddg(i,j,k,l)&T2s(i,j)(1,2,1,1)");
  test_for_zero(t4ddg_1(1,2,1,2) - (t4ddg_2(1,2,1,2)*t2s_3(1,2))
		,"T4ddg(i,j,k,l)&T2s(i,j)(1,2,1,2)");
  test_for_zero(t4ddg_1(1,2,2,0) - (t4ddg_2(1,2,2,0)*t2s_3(1,2))
		,"T4ddg(i,j,k,l)&T2s(i,j)(1,2,2,0)");
  test_for_zero(t4ddg_1(1,2,2,1) - (t4ddg_2(1,2,2,1)*t2s_3(1,2))
		,"T4ddg(i,j,k,l)&T2s(i,j)(1,2,2,1)");
  test_for_zero(t4ddg_1(1,2,2,2) - (t4ddg_2(1,2,2,2)*t2s_3(1,2))
		,"T4ddg(i,j,k,l)&T2s(i,j)(1,2,2,2)");
  test_for_zero(t4ddg_1(2,0,0,0) - (t4ddg_2(2,0,0,0)*t2s_3(2,0))
		,"T4ddg(i,j,k,l)&T2s(i,j)(2,0,0,0)");
  test_for_zero(t4ddg_1(2,0,0,1) - (t4ddg_2(2,0,0,1)*t2s_3(2,0))
		,"T4ddg(i,j,k,l)&T2s(i,j)(2,0,0,1)");
  test_for_zero(t4ddg_1(2,0,0,2) - (t4ddg_2(2,0,0,2)*t2s_3(2,0))
		,"T4ddg(i,j,k,l)&T2s(i,j)(2,0,0,2)");
  test_for_zero(t4ddg_1(2,0,1,0) - (t4ddg_2(2,0,1,0)*t2s_3(2,0))
		,"T4ddg(i,j,k,l)&T2s(i,j)(2,0,1,0)");
  test_for_zero(t4ddg_1(2,0,1,1) - (t4ddg_2(2,0,1,1)*t2s_3(2,0))
		,"T4ddg(i,j,k,l)&T2s(i,j)(2,0,1,1)");
  test_for_zero(t4ddg_1(2,0,1,2) - (t4ddg_2(2,0,1,2)*t2s_3(2,0))
		,"T4ddg(i,j,k,l)&T2s(i,j)(2,0,1,2)");
  test_for_zero(t4ddg_1(2,0,2,0) - (t4ddg_2(2,0,2,0)*t2s_3(2,0))
		,"T4ddg(i,j,k,l)&T2s(i,j)(2,0,2,0)");
  test_for_zero(t4ddg_1(2,0,2,1) - (t4ddg_2(2,0,2,1)*t2s_3(2,0))
		,"T4ddg(i,j,k,l)&T2s(i,j)(2,0,2,1)");
  test_for_zero(t4ddg_1(2,0,2,2) - (t4ddg_2(2,0,2,2)*t2s_3(2,0))
		,"T4ddg(i,j,k,l)&T2s(i,j)(2,0,2,2)");
  test_for_zero(t4ddg_1(2,1,0,0) - (t4ddg_2(2,1,0,0)*t2s_3(2,1))
		,"T4ddg(i,j,k,l)&T2s(i,j)(2,1,0,0)");
  test_for_zero(t4ddg_1(2,1,0,1) - (t4ddg_2(2,1,0,1)*t2s_3(2,1))
		,"T4ddg(i,j,k,l)&T2s(i,j)(2,1,0,1)");
  test_for_zero(t4ddg_1(2,1,0,2) - (t4ddg_2(2,1,0,2)*t2s_3(2,1))
		,"T4ddg(i,j,k,l)&T2s(i,j)(2,1,0,2)");
  test_for_zero(t4ddg_1(2,1,1,0) - (t4ddg_2(2,1,1,0)*t2s_3(2,1))
		,"T4ddg(i,j,k,l)&T2s(i,j)(2,1,1,0)");
  test_for_zero(t4ddg_1(2,1,1,1) - (t4ddg_2(2,1,1,1)*t2s_3(2,1))
		,"T4ddg(i,j,k,l)&T2s(i,j)(2,1,1,1)");
  test_for_zero(t4ddg_1(2,1,1,2) - (t4ddg_2(2,1,1,2)*t2s_3(2,1))
		,"T4ddg(i,j,k,l)&T2s(i,j)(2,1,1,2)");
  test_for_zero(t4ddg_1(2,1,2,0) - (t4ddg_2(2,1,2,0)*t2s_3(2,1))
		,"T4ddg(i,j,k,l)&T2s(i,j)(2,1,2,0)");
  test_for_zero(t4ddg_1(2,1,2,1) - (t4ddg_2(2,1,2,1)*t2s_3(2,1))
		,"T4ddg(i,j,k,l)&T2s(i,j)(2,1,2,1)");
  test_for_zero(t4ddg_1(2,1,2,2) - (t4ddg_2(2,1,2,2)*t2s_3(2,1))
		,"T4ddg(i,j,k,l)&T2s(i,j)(2,1,2,2)");
  test_for_zero(t4ddg_1(2,2,0,0) - (t4ddg_2(2,2,0,0)*t2s_3(2,2))
		,"T4ddg(i,j,k,l)&T2s(i,j)(2,2,0,0)");
  test_for_zero(t4ddg_1(2,2,0,1) - (t4ddg_2(2,2,0,1)*t2s_3(2,2))
		,"T4ddg(i,j,k,l)&T2s(i,j)(2,2,0,1)");
  test_for_zero(t4ddg_1(2,2,0,2) - (t4ddg_2(2,2,0,2)*t2s_3(2,2))
		,"T4ddg(i,j,k,l)&T2s(i,j)(2,2,0,2)");
  test_for_zero(t4ddg_1(2,2,1,0) - (t4ddg_2(2,2,1,0)*t2s_3(2,2))
		,"T4ddg(i,j,k,l)&T2s(i,j)(2,2,1,0)");
  test_for_zero(t4ddg_1(2,2,1,1) - (t4ddg_2(2,2,1,1)*t2s_3(2,2))
		,"T4ddg(i,j,k,l)&T2s(i,j)(2,2,1,1)");
  test_for_zero(t4ddg_1(2,2,1,2) - (t4ddg_2(2,2,1,2)*t2s_3(2,2))
		,"T4ddg(i,j,k,l)&T2s(i,j)(2,2,1,2)");
  test_for_zero(t4ddg_1(2,2,2,0) - (t4ddg_2(2,2,2,0)*t2s_3(2,2))
		,"T4ddg(i,j,k,l)&T2s(i,j)(2,2,2,0)");
  test_for_zero(t4ddg_1(2,2,2,1) - (t4ddg_2(2,2,2,1)*t2s_3(2,2))
		,"T4ddg(i,j,k,l)&T2s(i,j)(2,2,2,1)");
  test_for_zero(t4ddg_1(2,2,2,2) - (t4ddg_2(2,2,2,2)*t2s_3(2,2))
		,"T4ddg(i,j,k,l)&T2s(i,j)(2,2,2,2)");

  t4ddg_1(i,j,k,l)=(t2s_2(i,j)&t4ddg_2(i,j,k,l));
  test_for_zero(t4ddg_1(0,0,0,0) - (t4ddg_2(0,0,0,0)*t2s_2(0,0))
		,"T2s(i,j)&T4ddg(i,j,k,l)(0,0,0,0)");
  test_for_zero(t4ddg_1(0,0,0,1) - (t4ddg_2(0,0,0,1)*t2s_2(0,0))
		,"T2s(i,j)&T4ddg(i,j,k,l)(0,0,0,1)");
  test_for_zero(t4ddg_1(0,0,0,2) - (t4ddg_2(0,0,0,2)*t2s_2(0,0))
		,"T2s(i,j)&T4ddg(i,j,k,l)(0,0,0,2)");
  test_for_zero(t4ddg_1(0,0,1,0) - (t4ddg_2(0,0,1,0)*t2s_2(0,0))
		,"T2s(i,j)&T4ddg(i,j,k,l)(0,0,1,0)");
  test_for_zero(t4ddg_1(0,0,1,1) - (t4ddg_2(0,0,1,1)*t2s_2(0,0))
		,"T2s(i,j)&T4ddg(i,j,k,l)(0,0,1,1)");
  test_for_zero(t4ddg_1(0,0,1,2) - (t4ddg_2(0,0,1,2)*t2s_2(0,0))
		,"T2s(i,j)&T4ddg(i,j,k,l)(0,0,1,2)");
  test_for_zero(t4ddg_1(0,0,2,0) - (t4ddg_2(0,0,2,0)*t2s_2(0,0))
		,"T2s(i,j)&T4ddg(i,j,k,l)(0,0,2,0)");
  test_for_zero(t4ddg_1(0,0,2,1) - (t4ddg_2(0,0,2,1)*t2s_2(0,0))
		,"T2s(i,j)&T4ddg(i,j,k,l)(0,0,2,1)");
  test_for_zero(t4ddg_1(0,0,2,2) - (t4ddg_2(0,0,2,2)*t2s_2(0,0))
		,"T2s(i,j)&T4ddg(i,j,k,l)(0,0,2,2)");
  test_for_zero(t4ddg_1(0,1,0,0) - (t4ddg_2(0,1,0,0)*t2s_2(0,1))
		,"T2s(i,j)&T4ddg(i,j,k,l)(0,1,0,0)");
  test_for_zero(t4ddg_1(0,1,0,1) - (t4ddg_2(0,1,0,1)*t2s_2(0,1))
		,"T2s(i,j)&T4ddg(i,j,k,l)(0,1,0,1)");
  test_for_zero(t4ddg_1(0,1,0,2) - (t4ddg_2(0,1,0,2)*t2s_2(0,1))
		,"T2s(i,j)&T4ddg(i,j,k,l)(0,1,0,2)");
  test_for_zero(t4ddg_1(0,1,1,0) - (t4ddg_2(0,1,1,0)*t2s_2(0,1))
		,"T2s(i,j)&T4ddg(i,j,k,l)(0,1,1,0)");
  test_for_zero(t4ddg_1(0,1,1,1) - (t4ddg_2(0,1,1,1)*t2s_2(0,1))
		,"T2s(i,j)&T4ddg(i,j,k,l)(0,1,1,1)");
  test_for_zero(t4ddg_1(0,1,1,2) - (t4ddg_2(0,1,1,2)*t2s_2(0,1))
		,"T2s(i,j)&T4ddg(i,j,k,l)(0,1,1,2)");
  test_for_zero(t4ddg_1(0,1,2,0) - (t4ddg_2(0,1,2,0)*t2s_2(0,1))
		,"T2s(i,j)&T4ddg(i,j,k,l)(0,1,2,0)");
  test_for_zero(t4ddg_1(0,1,2,1) - (t4ddg_2(0,1,2,1)*t2s_2(0,1))
		,"T2s(i,j)&T4ddg(i,j,k,l)(0,1,2,1)");
  test_for_zero(t4ddg_1(0,1,2,2) - (t4ddg_2(0,1,2,2)*t2s_2(0,1))
		,"T2s(i,j)&T4ddg(i,j,k,l)(0,1,2,2)");
  test_for_zero(t4ddg_1(0,2,0,0) - (t4ddg_2(0,2,0,0)*t2s_2(0,2))
		,"T2s(i,j)&T4ddg(i,j,k,l)(0,2,0,0)");
  test_for_zero(t4ddg_1(0,2,0,1) - (t4ddg_2(0,2,0,1)*t2s_2(0,2))
		,"T2s(i,j)&T4ddg(i,j,k,l)(0,2,0,1)");
  test_for_zero(t4ddg_1(0,2,0,2) - (t4ddg_2(0,2,0,2)*t2s_2(0,2))
		,"T2s(i,j)&T4ddg(i,j,k,l)(0,2,0,2)");
  test_for_zero(t4ddg_1(0,2,1,0) - (t4ddg_2(0,2,1,0)*t2s_2(0,2))
		,"T2s(i,j)&T4ddg(i,j,k,l)(0,2,1,0)");
  test_for_zero(t4ddg_1(0,2,1,1) - (t4ddg_2(0,2,1,1)*t2s_2(0,2))
		,"T2s(i,j)&T4ddg(i,j,k,l)(0,2,1,1)");
  test_for_zero(t4ddg_1(0,2,1,2) - (t4ddg_2(0,2,1,2)*t2s_2(0,2))
		,"T2s(i,j)&T4ddg(i,j,k,l)(0,2,1,2)");
  test_for_zero(t4ddg_1(0,2,2,0) - (t4ddg_2(0,2,2,0)*t2s_2(0,2))
		,"T2s(i,j)&T4ddg(i,j,k,l)(0,2,2,0)");
  test_for_zero(t4ddg_1(0,2,2,1) - (t4ddg_2(0,2,2,1)*t2s_2(0,2))
		,"T2s(i,j)&T4ddg(i,j,k,l)(0,2,2,1)");
  test_for_zero(t4ddg_1(0,2,2,2) - (t4ddg_2(0,2,2,2)*t2s_2(0,2))
		,"T2s(i,j)&T4ddg(i,j,k,l)(0,2,2,2)");
  test_for_zero(t4ddg_1(1,0,0,0) - (t4ddg_2(1,0,0,0)*t2s_2(1,0))
		,"T2s(i,j)&T4ddg(i,j,k,l)(1,0,0,0)");
  test_for_zero(t4ddg_1(1,0,0,1) - (t4ddg_2(1,0,0,1)*t2s_2(1,0))
		,"T2s(i,j)&T4ddg(i,j,k,l)(1,0,0,1)");
  test_for_zero(t4ddg_1(1,0,0,2) - (t4ddg_2(1,0,0,2)*t2s_2(1,0))
		,"T2s(i,j)&T4ddg(i,j,k,l)(1,0,0,2)");
  test_for_zero(t4ddg_1(1,0,1,0) - (t4ddg_2(1,0,1,0)*t2s_2(1,0))
		,"T2s(i,j)&T4ddg(i,j,k,l)(1,0,1,0)");
  test_for_zero(t4ddg_1(1,0,1,1) - (t4ddg_2(1,0,1,1)*t2s_2(1,0))
		,"T2s(i,j)&T4ddg(i,j,k,l)(1,0,1,1)");
  test_for_zero(t4ddg_1(1,0,1,2) - (t4ddg_2(1,0,1,2)*t2s_2(1,0))
		,"T2s(i,j)&T4ddg(i,j,k,l)(1,0,1,2)");
  test_for_zero(t4ddg_1(1,0,2,0) - (t4ddg_2(1,0,2,0)*t2s_2(1,0))
		,"T2s(i,j)&T4ddg(i,j,k,l)(1,0,2,0)");
  test_for_zero(t4ddg_1(1,0,2,1) - (t4ddg_2(1,0,2,1)*t2s_2(1,0))
		,"T2s(i,j)&T4ddg(i,j,k,l)(1,0,2,1)");
  test_for_zero(t4ddg_1(1,0,2,2) - (t4ddg_2(1,0,2,2)*t2s_2(1,0))
		,"T2s(i,j)&T4ddg(i,j,k,l)(1,0,2,2)");
  test_for_zero(t4ddg_1(1,1,0,0) - (t4ddg_2(1,1,0,0)*t2s_2(1,1))
		,"T2s(i,j)&T4ddg(i,j,k,l)(1,1,0,0)");
  test_for_zero(t4ddg_1(1,1,0,1) - (t4ddg_2(1,1,0,1)*t2s_2(1,1))
		,"T2s(i,j)&T4ddg(i,j,k,l)(1,1,0,1)");
  test_for_zero(t4ddg_1(1,1,0,2) - (t4ddg_2(1,1,0,2)*t2s_2(1,1))
		,"T2s(i,j)&T4ddg(i,j,k,l)(1,1,0,2)");
  test_for_zero(t4ddg_1(1,1,1,0) - (t4ddg_2(1,1,1,0)*t2s_2(1,1))
		,"T2s(i,j)&T4ddg(i,j,k,l)(1,1,1,0)");
  test_for_zero(t4ddg_1(1,1,1,1) - (t4ddg_2(1,1,1,1)*t2s_2(1,1))
		,"T2s(i,j)&T4ddg(i,j,k,l)(1,1,1,1)");
  test_for_zero(t4ddg_1(1,1,1,2) - (t4ddg_2(1,1,1,2)*t2s_2(1,1))
		,"T2s(i,j)&T4ddg(i,j,k,l)(1,1,1,2)");
  test_for_zero(t4ddg_1(1,1,2,0) - (t4ddg_2(1,1,2,0)*t2s_2(1,1))
		,"T2s(i,j)&T4ddg(i,j,k,l)(1,1,2,0)");
  test_for_zero(t4ddg_1(1,1,2,1) - (t4ddg_2(1,1,2,1)*t2s_2(1,1))
		,"T2s(i,j)&T4ddg(i,j,k,l)(1,1,2,1)");
  test_for_zero(t4ddg_1(1,1,2,2) - (t4ddg_2(1,1,2,2)*t2s_2(1,1))
		,"T2s(i,j)&T4ddg(i,j,k,l)(1,1,2,2)");
  test_for_zero(t4ddg_1(1,2,0,0) - (t4ddg_2(1,2,0,0)*t2s_2(1,2))
		,"T2s(i,j)&T4ddg(i,j,k,l)(1,2,0,0)");
  test_for_zero(t4ddg_1(1,2,0,1) - (t4ddg_2(1,2,0,1)*t2s_2(1,2))
		,"T2s(i,j)&T4ddg(i,j,k,l)(1,2,0,1)");
  test_for_zero(t4ddg_1(1,2,0,2) - (t4ddg_2(1,2,0,2)*t2s_2(1,2))
		,"T2s(i,j)&T4ddg(i,j,k,l)(1,2,0,2)");
  test_for_zero(t4ddg_1(1,2,1,0) - (t4ddg_2(1,2,1,0)*t2s_2(1,2))
		,"T2s(i,j)&T4ddg(i,j,k,l)(1,2,1,0)");
  test_for_zero(t4ddg_1(1,2,1,1) - (t4ddg_2(1,2,1,1)*t2s_2(1,2))
		,"T2s(i,j)&T4ddg(i,j,k,l)(1,2,1,1)");
  test_for_zero(t4ddg_1(1,2,1,2) - (t4ddg_2(1,2,1,2)*t2s_2(1,2))
		,"T2s(i,j)&T4ddg(i,j,k,l)(1,2,1,2)");
  test_for_zero(t4ddg_1(1,2,2,0) - (t4ddg_2(1,2,2,0)*t2s_2(1,2))
		,"T2s(i,j)&T4ddg(i,j,k,l)(1,2,2,0)");
  test_for_zero(t4ddg_1(1,2,2,1) - (t4ddg_2(1,2,2,1)*t2s_2(1,2))
		,"T2s(i,j)&T4ddg(i,j,k,l)(1,2,2,1)");
  test_for_zero(t4ddg_1(1,2,2,2) - (t4ddg_2(1,2,2,2)*t2s_2(1,2))
		,"T2s(i,j)&T4ddg(i,j,k,l)(1,2,2,2)");
  test_for_zero(t4ddg_1(2,0,0,0) - (t4ddg_2(2,0,0,0)*t2s_2(2,0))
		,"T2s(i,j)&T4ddg(i,j,k,l)(2,0,0,0)");
  test_for_zero(t4ddg_1(2,0,0,1) - (t4ddg_2(2,0,0,1)*t2s_2(2,0))
		,"T2s(i,j)&T4ddg(i,j,k,l)(2,0,0,1)");
  test_for_zero(t4ddg_1(2,0,0,2) - (t4ddg_2(2,0,0,2)*t2s_2(2,0))
		,"T2s(i,j)&T4ddg(i,j,k,l)(2,0,0,2)");
  test_for_zero(t4ddg_1(2,0,1,0) - (t4ddg_2(2,0,1,0)*t2s_2(2,0))
		,"T2s(i,j)&T4ddg(i,j,k,l)(2,0,1,0)");
  test_for_zero(t4ddg_1(2,0,1,1) - (t4ddg_2(2,0,1,1)*t2s_2(2,0))
		,"T2s(i,j)&T4ddg(i,j,k,l)(2,0,1,1)");
  test_for_zero(t4ddg_1(2,0,1,2) - (t4ddg_2(2,0,1,2)*t2s_2(2,0))
		,"T2s(i,j)&T4ddg(i,j,k,l)(2,0,1,2)");
  test_for_zero(t4ddg_1(2,0,2,0) - (t4ddg_2(2,0,2,0)*t2s_2(2,0))
		,"T2s(i,j)&T4ddg(i,j,k,l)(2,0,2,0)");
  test_for_zero(t4ddg_1(2,0,2,1) - (t4ddg_2(2,0,2,1)*t2s_2(2,0))
		,"T2s(i,j)&T4ddg(i,j,k,l)(2,0,2,1)");
  test_for_zero(t4ddg_1(2,0,2,2) - (t4ddg_2(2,0,2,2)*t2s_2(2,0))
		,"T2s(i,j)&T4ddg(i,j,k,l)(2,0,2,2)");
  test_for_zero(t4ddg_1(2,1,0,0) - (t4ddg_2(2,1,0,0)*t2s_2(2,1))
		,"T2s(i,j)&T4ddg(i,j,k,l)(2,1,0,0)");
  test_for_zero(t4ddg_1(2,1,0,1) - (t4ddg_2(2,1,0,1)*t2s_2(2,1))
		,"T2s(i,j)&T4ddg(i,j,k,l)(2,1,0,1)");
  test_for_zero(t4ddg_1(2,1,0,2) - (t4ddg_2(2,1,0,2)*t2s_2(2,1))
		,"T2s(i,j)&T4ddg(i,j,k,l)(2,1,0,2)");
  test_for_zero(t4ddg_1(2,1,1,0) - (t4ddg_2(2,1,1,0)*t2s_2(2,1))
		,"T2s(i,j)&T4ddg(i,j,k,l)(2,1,1,0)");
  test_for_zero(t4ddg_1(2,1,1,1) - (t4ddg_2(2,1,1,1)*t2s_2(2,1))
		,"T2s(i,j)&T4ddg(i,j,k,l)(2,1,1,1)");
  test_for_zero(t4ddg_1(2,1,1,2) - (t4ddg_2(2,1,1,2)*t2s_2(2,1))
		,"T2s(i,j)&T4ddg(i,j,k,l)(2,1,1,2)");
  test_for_zero(t4ddg_1(2,1,2,0) - (t4ddg_2(2,1,2,0)*t2s_2(2,1))
		,"T2s(i,j)&T4ddg(i,j,k,l)(2,1,2,0)");
  test_for_zero(t4ddg_1(2,1,2,1) - (t4ddg_2(2,1,2,1)*t2s_2(2,1))
		,"T2s(i,j)&T4ddg(i,j,k,l)(2,1,2,1)");
  test_for_zero(t4ddg_1(2,1,2,2) - (t4ddg_2(2,1,2,2)*t2s_2(2,1))
		,"T2s(i,j)&T4ddg(i,j,k,l)(2,1,2,2)");
  test_for_zero(t4ddg_1(2,2,0,0) - (t4ddg_2(2,2,0,0)*t2s_2(2,2))
		,"T2s(i,j)&T4ddg(i,j,k,l)(2,2,0,0)");
  test_for_zero(t4ddg_1(2,2,0,1) - (t4ddg_2(2,2,0,1)*t2s_2(2,2))
		,"T2s(i,j)&T4ddg(i,j,k,l)(2,2,0,1)");
  test_for_zero(t4ddg_1(2,2,0,2) - (t4ddg_2(2,2,0,2)*t2s_2(2,2))
		,"T2s(i,j)&T4ddg(i,j,k,l)(2,2,0,2)");
  test_for_zero(t4ddg_1(2,2,1,0) - (t4ddg_2(2,2,1,0)*t2s_2(2,2))
		,"T2s(i,j)&T4ddg(i,j,k,l)(2,2,1,0)");
  test_for_zero(t4ddg_1(2,2,1,1) - (t4ddg_2(2,2,1,1)*t2s_2(2,2))
		,"T2s(i,j)&T4ddg(i,j,k,l)(2,2,1,1)");
  test_for_zero(t4ddg_1(2,2,1,2) - (t4ddg_2(2,2,1,2)*t2s_2(2,2))
		,"T2s(i,j)&T4ddg(i,j,k,l)(2,2,1,2)");
  test_for_zero(t4ddg_1(2,2,2,0) - (t4ddg_2(2,2,2,0)*t2s_2(2,2))
		,"T2s(i,j)&T4ddg(i,j,k,l)(2,2,2,0)");
  test_for_zero(t4ddg_1(2,2,2,1) - (t4ddg_2(2,2,2,1)*t2s_2(2,2))
		,"T2s(i,j)&T4ddg(i,j,k,l)(2,2,2,1)");
  test_for_zero(t4ddg_1(2,2,2,2) - (t4ddg_2(2,2,2,2)*t2s_2(2,2))
		,"T2s(i,j)&T4ddg(i,j,k,l)(2,2,2,2)");

  t4ddg_1(i,j,k,l)=(t4ddg_2(i,j,k,l)&t2s_3(k,l));
  test_for_zero(t4ddg_1(0,0,0,0) - (t4ddg_2(0,0,0,0)*t2s_3(0,0))
		,"T4ddg(i,j,k,l)&T2s(k,l)(0,0,0,0)");
  test_for_zero(t4ddg_1(0,0,0,1) - (t4ddg_2(0,0,0,1)*t2s_3(0,1))
		,"T4ddg(i,j,k,l)&T2s(k,l)(0,0,0,1)");
  test_for_zero(t4ddg_1(0,0,0,2) - (t4ddg_2(0,0,0,2)*t2s_3(0,2))
		,"T4ddg(i,j,k,l)&T2s(k,l)(0,0,0,2)");
  test_for_zero(t4ddg_1(0,0,1,0) - (t4ddg_2(0,0,1,0)*t2s_3(1,0))
		,"T4ddg(i,j,k,l)&T2s(k,l)(0,0,1,0)");
  test_for_zero(t4ddg_1(0,0,1,1) - (t4ddg_2(0,0,1,1)*t2s_3(1,1))
		,"T4ddg(i,j,k,l)&T2s(k,l)(0,0,1,1)");
  test_for_zero(t4ddg_1(0,0,1,2) - (t4ddg_2(0,0,1,2)*t2s_3(1,2))
		,"T4ddg(i,j,k,l)&T2s(k,l)(0,0,1,2)");
  test_for_zero(t4ddg_1(0,0,2,0) - (t4ddg_2(0,0,2,0)*t2s_3(2,0))
		,"T4ddg(i,j,k,l)&T2s(k,l)(0,0,2,0)");
  test_for_zero(t4ddg_1(0,0,2,1) - (t4ddg_2(0,0,2,1)*t2s_3(2,1))
		,"T4ddg(i,j,k,l)&T2s(k,l)(0,0,2,1)");
  test_for_zero(t4ddg_1(0,0,2,2) - (t4ddg_2(0,0,2,2)*t2s_3(2,2))
		,"T4ddg(i,j,k,l)&T2s(k,l)(0,0,2,2)");
  test_for_zero(t4ddg_1(0,1,0,0) - (t4ddg_2(0,1,0,0)*t2s_3(0,0))
		,"T4ddg(i,j,k,l)&T2s(k,l)(0,1,0,0)");
  test_for_zero(t4ddg_1(0,1,0,1) - (t4ddg_2(0,1,0,1)*t2s_3(0,1))
		,"T4ddg(i,j,k,l)&T2s(k,l)(0,1,0,1)");
  test_for_zero(t4ddg_1(0,1,0,2) - (t4ddg_2(0,1,0,2)*t2s_3(0,2))
		,"T4ddg(i,j,k,l)&T2s(k,l)(0,1,0,2)");
  test_for_zero(t4ddg_1(0,1,1,0) - (t4ddg_2(0,1,1,0)*t2s_3(1,0))
		,"T4ddg(i,j,k,l)&T2s(k,l)(0,1,1,0)");
  test_for_zero(t4ddg_1(0,1,1,1) - (t4ddg_2(0,1,1,1)*t2s_3(1,1))
		,"T4ddg(i,j,k,l)&T2s(k,l)(0,1,1,1)");
  test_for_zero(t4ddg_1(0,1,1,2) - (t4ddg_2(0,1,1,2)*t2s_3(1,2))
		,"T4ddg(i,j,k,l)&T2s(k,l)(0,1,1,2)");
  test_for_zero(t4ddg_1(0,1,2,0) - (t4ddg_2(0,1,2,0)*t2s_3(2,0))
		,"T4ddg(i,j,k,l)&T2s(k,l)(0,1,2,0)");
  test_for_zero(t4ddg_1(0,1,2,1) - (t4ddg_2(0,1,2,1)*t2s_3(2,1))
		,"T4ddg(i,j,k,l)&T2s(k,l)(0,1,2,1)");
  test_for_zero(t4ddg_1(0,1,2,2) - (t4ddg_2(0,1,2,2)*t2s_3(2,2))
		,"T4ddg(i,j,k,l)&T2s(k,l)(0,1,2,2)");
  test_for_zero(t4ddg_1(0,2,0,0) - (t4ddg_2(0,2,0,0)*t2s_3(0,0))
		,"T4ddg(i,j,k,l)&T2s(k,l)(0,2,0,0)");
  test_for_zero(t4ddg_1(0,2,0,1) - (t4ddg_2(0,2,0,1)*t2s_3(0,1))
		,"T4ddg(i,j,k,l)&T2s(k,l)(0,2,0,1)");
  test_for_zero(t4ddg_1(0,2,0,2) - (t4ddg_2(0,2,0,2)*t2s_3(0,2))
		,"T4ddg(i,j,k,l)&T2s(k,l)(0,2,0,2)");
  test_for_zero(t4ddg_1(0,2,1,0) - (t4ddg_2(0,2,1,0)*t2s_3(1,0))
		,"T4ddg(i,j,k,l)&T2s(k,l)(0,2,1,0)");
  test_for_zero(t4ddg_1(0,2,1,1) - (t4ddg_2(0,2,1,1)*t2s_3(1,1))
		,"T4ddg(i,j,k,l)&T2s(k,l)(0,2,1,1)");
  test_for_zero(t4ddg_1(0,2,1,2) - (t4ddg_2(0,2,1,2)*t2s_3(1,2))
		,"T4ddg(i,j,k,l)&T2s(k,l)(0,2,1,2)");
  test_for_zero(t4ddg_1(0,2,2,0) - (t4ddg_2(0,2,2,0)*t2s_3(2,0))
		,"T4ddg(i,j,k,l)&T2s(k,l)(0,2,2,0)");
  test_for_zero(t4ddg_1(0,2,2,1) - (t4ddg_2(0,2,2,1)*t2s_3(2,1))
		,"T4ddg(i,j,k,l)&T2s(k,l)(0,2,2,1)");
  test_for_zero(t4ddg_1(0,2,2,2) - (t4ddg_2(0,2,2,2)*t2s_3(2,2))
		,"T4ddg(i,j,k,l)&T2s(k,l)(0,2,2,2)");
  test_for_zero(t4ddg_1(1,0,0,0) - (t4ddg_2(1,0,0,0)*t2s_3(0,0))
		,"T4ddg(i,j,k,l)&T2s(k,l)(1,0,0,0)");
  test_for_zero(t4ddg_1(1,0,0,1) - (t4ddg_2(1,0,0,1)*t2s_3(0,1))
		,"T4ddg(i,j,k,l)&T2s(k,l)(1,0,0,1)");
  test_for_zero(t4ddg_1(1,0,0,2) - (t4ddg_2(1,0,0,2)*t2s_3(0,2))
		,"T4ddg(i,j,k,l)&T2s(k,l)(1,0,0,2)");
  test_for_zero(t4ddg_1(1,0,1,0) - (t4ddg_2(1,0,1,0)*t2s_3(1,0))
		,"T4ddg(i,j,k,l)&T2s(k,l)(1,0,1,0)");
  test_for_zero(t4ddg_1(1,0,1,1) - (t4ddg_2(1,0,1,1)*t2s_3(1,1))
		,"T4ddg(i,j,k,l)&T2s(k,l)(1,0,1,1)");
  test_for_zero(t4ddg_1(1,0,1,2) - (t4ddg_2(1,0,1,2)*t2s_3(1,2))
		,"T4ddg(i,j,k,l)&T2s(k,l)(1,0,1,2)");
  test_for_zero(t4ddg_1(1,0,2,0) - (t4ddg_2(1,0,2,0)*t2s_3(2,0))
		,"T4ddg(i,j,k,l)&T2s(k,l)(1,0,2,0)");
  test_for_zero(t4ddg_1(1,0,2,1) - (t4ddg_2(1,0,2,1)*t2s_3(2,1))
		,"T4ddg(i,j,k,l)&T2s(k,l)(1,0,2,1)");
  test_for_zero(t4ddg_1(1,0,2,2) - (t4ddg_2(1,0,2,2)*t2s_3(2,2))
		,"T4ddg(i,j,k,l)&T2s(k,l)(1,0,2,2)");
  test_for_zero(t4ddg_1(1,1,0,0) - (t4ddg_2(1,1,0,0)*t2s_3(0,0))
		,"T4ddg(i,j,k,l)&T2s(k,l)(1,1,0,0)");
  test_for_zero(t4ddg_1(1,1,0,1) - (t4ddg_2(1,1,0,1)*t2s_3(0,1))
		,"T4ddg(i,j,k,l)&T2s(k,l)(1,1,0,1)");
  test_for_zero(t4ddg_1(1,1,0,2) - (t4ddg_2(1,1,0,2)*t2s_3(0,2))
		,"T4ddg(i,j,k,l)&T2s(k,l)(1,1,0,2)");
  test_for_zero(t4ddg_1(1,1,1,0) - (t4ddg_2(1,1,1,0)*t2s_3(1,0))
		,"T4ddg(i,j,k,l)&T2s(k,l)(1,1,1,0)");
  test_for_zero(t4ddg_1(1,1,1,1) - (t4ddg_2(1,1,1,1)*t2s_3(1,1))
		,"T4ddg(i,j,k,l)&T2s(k,l)(1,1,1,1)");
  test_for_zero(t4ddg_1(1,1,1,2) - (t4ddg_2(1,1,1,2)*t2s_3(1,2))
		,"T4ddg(i,j,k,l)&T2s(k,l)(1,1,1,2)");
  test_for_zero(t4ddg_1(1,1,2,0) - (t4ddg_2(1,1,2,0)*t2s_3(2,0))
		,"T4ddg(i,j,k,l)&T2s(k,l)(1,1,2,0)");
  test_for_zero(t4ddg_1(1,1,2,1) - (t4ddg_2(1,1,2,1)*t2s_3(2,1))
		,"T4ddg(i,j,k,l)&T2s(k,l)(1,1,2,1)");
  test_for_zero(t4ddg_1(1,1,2,2) - (t4ddg_2(1,1,2,2)*t2s_3(2,2))
		,"T4ddg(i,j,k,l)&T2s(k,l)(1,1,2,2)");
  test_for_zero(t4ddg_1(1,2,0,0) - (t4ddg_2(1,2,0,0)*t2s_3(0,0))
		,"T4ddg(i,j,k,l)&T2s(k,l)(1,2,0,0)");
  test_for_zero(t4ddg_1(1,2,0,1) - (t4ddg_2(1,2,0,1)*t2s_3(0,1))
		,"T4ddg(i,j,k,l)&T2s(k,l)(1,2,0,1)");
  test_for_zero(t4ddg_1(1,2,0,2) - (t4ddg_2(1,2,0,2)*t2s_3(0,2))
		,"T4ddg(i,j,k,l)&T2s(k,l)(1,2,0,2)");
  test_for_zero(t4ddg_1(1,2,1,0) - (t4ddg_2(1,2,1,0)*t2s_3(1,0))
		,"T4ddg(i,j,k,l)&T2s(k,l)(1,2,1,0)");
  test_for_zero(t4ddg_1(1,2,1,1) - (t4ddg_2(1,2,1,1)*t2s_3(1,1))
		,"T4ddg(i,j,k,l)&T2s(k,l)(1,2,1,1)");
  test_for_zero(t4ddg_1(1,2,1,2) - (t4ddg_2(1,2,1,2)*t2s_3(1,2))
		,"T4ddg(i,j,k,l)&T2s(k,l)(1,2,1,2)");
  test_for_zero(t4ddg_1(1,2,2,0) - (t4ddg_2(1,2,2,0)*t2s_3(2,0))
		,"T4ddg(i,j,k,l)&T2s(k,l)(1,2,2,0)");
  test_for_zero(t4ddg_1(1,2,2,1) - (t4ddg_2(1,2,2,1)*t2s_3(2,1))
		,"T4ddg(i,j,k,l)&T2s(k,l)(1,2,2,1)");
  test_for_zero(t4ddg_1(1,2,2,2) - (t4ddg_2(1,2,2,2)*t2s_3(2,2))
		,"T4ddg(i,j,k,l)&T2s(k,l)(1,2,2,2)");
  test_for_zero(t4ddg_1(2,0,0,0) - (t4ddg_2(2,0,0,0)*t2s_3(0,0))
		,"T4ddg(i,j,k,l)&T2s(k,l)(2,0,0,0)");
  test_for_zero(t4ddg_1(2,0,0,1) - (t4ddg_2(2,0,0,1)*t2s_3(0,1))
		,"T4ddg(i,j,k,l)&T2s(k,l)(2,0,0,1)");
  test_for_zero(t4ddg_1(2,0,0,2) - (t4ddg_2(2,0,0,2)*t2s_3(0,2))
		,"T4ddg(i,j,k,l)&T2s(k,l)(2,0,0,2)");
  test_for_zero(t4ddg_1(2,0,1,0) - (t4ddg_2(2,0,1,0)*t2s_3(1,0))
		,"T4ddg(i,j,k,l)&T2s(k,l)(2,0,1,0)");
  test_for_zero(t4ddg_1(2,0,1,1) - (t4ddg_2(2,0,1,1)*t2s_3(1,1))
		,"T4ddg(i,j,k,l)&T2s(k,l)(2,0,1,1)");
  test_for_zero(t4ddg_1(2,0,1,2) - (t4ddg_2(2,0,1,2)*t2s_3(1,2))
		,"T4ddg(i,j,k,l)&T2s(k,l)(2,0,1,2)");
  test_for_zero(t4ddg_1(2,0,2,0) - (t4ddg_2(2,0,2,0)*t2s_3(2,0))
		,"T4ddg(i,j,k,l)&T2s(k,l)(2,0,2,0)");
  test_for_zero(t4ddg_1(2,0,2,1) - (t4ddg_2(2,0,2,1)*t2s_3(2,1))
		,"T4ddg(i,j,k,l)&T2s(k,l)(2,0,2,1)");
  test_for_zero(t4ddg_1(2,0,2,2) - (t4ddg_2(2,0,2,2)*t2s_3(2,2))
		,"T4ddg(i,j,k,l)&T2s(k,l)(2,0,2,2)");
  test_for_zero(t4ddg_1(2,1,0,0) - (t4ddg_2(2,1,0,0)*t2s_3(0,0))
		,"T4ddg(i,j,k,l)&T2s(k,l)(2,1,0,0)");
  test_for_zero(t4ddg_1(2,1,0,1) - (t4ddg_2(2,1,0,1)*t2s_3(0,1))
		,"T4ddg(i,j,k,l)&T2s(k,l)(2,1,0,1)");
  test_for_zero(t4ddg_1(2,1,0,2) - (t4ddg_2(2,1,0,2)*t2s_3(0,2))
		,"T4ddg(i,j,k,l)&T2s(k,l)(2,1,0,2)");
  test_for_zero(t4ddg_1(2,1,1,0) - (t4ddg_2(2,1,1,0)*t2s_3(1,0))
		,"T4ddg(i,j,k,l)&T2s(k,l)(2,1,1,0)");
  test_for_zero(t4ddg_1(2,1,1,1) - (t4ddg_2(2,1,1,1)*t2s_3(1,1))
		,"T4ddg(i,j,k,l)&T2s(k,l)(2,1,1,1)");
  test_for_zero(t4ddg_1(2,1,1,2) - (t4ddg_2(2,1,1,2)*t2s_3(1,2))
		,"T4ddg(i,j,k,l)&T2s(k,l)(2,1,1,2)");
  test_for_zero(t4ddg_1(2,1,2,0) - (t4ddg_2(2,1,2,0)*t2s_3(2,0))
		,"T4ddg(i,j,k,l)&T2s(k,l)(2,1,2,0)");
  test_for_zero(t4ddg_1(2,1,2,1) - (t4ddg_2(2,1,2,1)*t2s_3(2,1))
		,"T4ddg(i,j,k,l)&T2s(k,l)(2,1,2,1)");
  test_for_zero(t4ddg_1(2,1,2,2) - (t4ddg_2(2,1,2,2)*t2s_3(2,2))
		,"T4ddg(i,j,k,l)&T2s(k,l)(2,1,2,2)");
  test_for_zero(t4ddg_1(2,2,0,0) - (t4ddg_2(2,2,0,0)*t2s_3(0,0))
		,"T4ddg(i,j,k,l)&T2s(k,l)(2,2,0,0)");
  test_for_zero(t4ddg_1(2,2,0,1) - (t4ddg_2(2,2,0,1)*t2s_3(0,1))
		,"T4ddg(i,j,k,l)&T2s(k,l)(2,2,0,1)");
  test_for_zero(t4ddg_1(2,2,0,2) - (t4ddg_2(2,2,0,2)*t2s_3(0,2))
		,"T4ddg(i,j,k,l)&T2s(k,l)(2,2,0,2)");
  test_for_zero(t4ddg_1(2,2,1,0) - (t4ddg_2(2,2,1,0)*t2s_3(1,0))
		,"T4ddg(i,j,k,l)&T2s(k,l)(2,2,1,0)");
  test_for_zero(t4ddg_1(2,2,1,1) - (t4ddg_2(2,2,1,1)*t2s_3(1,1))
		,"T4ddg(i,j,k,l)&T2s(k,l)(2,2,1,1)");
  test_for_zero(t4ddg_1(2,2,1,2) - (t4ddg_2(2,2,1,2)*t2s_3(1,2))
		,"T4ddg(i,j,k,l)&T2s(k,l)(2,2,1,2)");
  test_for_zero(t4ddg_1(2,2,2,0) - (t4ddg_2(2,2,2,0)*t2s_3(2,0))
		,"T4ddg(i,j,k,l)&T2s(k,l)(2,2,2,0)");
  test_for_zero(t4ddg_1(2,2,2,1) - (t4ddg_2(2,2,2,1)*t2s_3(2,1))
		,"T4ddg(i,j,k,l)&T2s(k,l)(2,2,2,1)");
  test_for_zero(t4ddg_1(2,2,2,2) - (t4ddg_2(2,2,2,2)*t2s_3(2,2))
		,"T4ddg(i,j,k,l)&T2s(k,l)(2,2,2,2)");

  t4ddg_1(i,j,k,l)=(t2s_2(k,l)&t4ddg_2(i,j,k,l));
  test_for_zero(t4ddg_1(0,0,0,0) - (t4ddg_2(0,0,0,0)*t2s_2(0,0))
		,"T2s(k,l)&T4ddg(i,j,k,l)(0,0,0,0)");
  test_for_zero(t4ddg_1(0,0,0,1) - (t4ddg_2(0,0,0,1)*t2s_2(0,1))
		,"T2s(k,l)&T4ddg(i,j,k,l)(0,0,0,1)");
  test_for_zero(t4ddg_1(0,0,0,2) - (t4ddg_2(0,0,0,2)*t2s_2(0,2))
		,"T2s(k,l)&T4ddg(i,j,k,l)(0,0,0,2)");
  test_for_zero(t4ddg_1(0,0,1,0) - (t4ddg_2(0,0,1,0)*t2s_2(1,0))
		,"T2s(k,l)&T4ddg(i,j,k,l)(0,0,1,0)");
  test_for_zero(t4ddg_1(0,0,1,1) - (t4ddg_2(0,0,1,1)*t2s_2(1,1))
		,"T2s(k,l)&T4ddg(i,j,k,l)(0,0,1,1)");
  test_for_zero(t4ddg_1(0,0,1,2) - (t4ddg_2(0,0,1,2)*t2s_2(1,2))
		,"T2s(k,l)&T4ddg(i,j,k,l)(0,0,1,2)");
  test_for_zero(t4ddg_1(0,0,2,0) - (t4ddg_2(0,0,2,0)*t2s_2(2,0))
		,"T2s(k,l)&T4ddg(i,j,k,l)(0,0,2,0)");
  test_for_zero(t4ddg_1(0,0,2,1) - (t4ddg_2(0,0,2,1)*t2s_2(2,1))
		,"T2s(k,l)&T4ddg(i,j,k,l)(0,0,2,1)");
  test_for_zero(t4ddg_1(0,0,2,2) - (t4ddg_2(0,0,2,2)*t2s_2(2,2))
		,"T2s(k,l)&T4ddg(i,j,k,l)(0,0,2,2)");
  test_for_zero(t4ddg_1(0,1,0,0) - (t4ddg_2(0,1,0,0)*t2s_2(0,0))
		,"T2s(k,l)&T4ddg(i,j,k,l)(0,1,0,0)");
  test_for_zero(t4ddg_1(0,1,0,1) - (t4ddg_2(0,1,0,1)*t2s_2(0,1))
		,"T2s(k,l)&T4ddg(i,j,k,l)(0,1,0,1)");
  test_for_zero(t4ddg_1(0,1,0,2) - (t4ddg_2(0,1,0,2)*t2s_2(0,2))
		,"T2s(k,l)&T4ddg(i,j,k,l)(0,1,0,2)");
  test_for_zero(t4ddg_1(0,1,1,0) - (t4ddg_2(0,1,1,0)*t2s_2(1,0))
		,"T2s(k,l)&T4ddg(i,j,k,l)(0,1,1,0)");
  test_for_zero(t4ddg_1(0,1,1,1) - (t4ddg_2(0,1,1,1)*t2s_2(1,1))
		,"T2s(k,l)&T4ddg(i,j,k,l)(0,1,1,1)");
  test_for_zero(t4ddg_1(0,1,1,2) - (t4ddg_2(0,1,1,2)*t2s_2(1,2))
		,"T2s(k,l)&T4ddg(i,j,k,l)(0,1,1,2)");
  test_for_zero(t4ddg_1(0,1,2,0) - (t4ddg_2(0,1,2,0)*t2s_2(2,0))
		,"T2s(k,l)&T4ddg(i,j,k,l)(0,1,2,0)");
  test_for_zero(t4ddg_1(0,1,2,1) - (t4ddg_2(0,1,2,1)*t2s_2(2,1))
		,"T2s(k,l)&T4ddg(i,j,k,l)(0,1,2,1)");
  test_for_zero(t4ddg_1(0,1,2,2) - (t4ddg_2(0,1,2,2)*t2s_2(2,2))
		,"T2s(k,l)&T4ddg(i,j,k,l)(0,1,2,2)");
  test_for_zero(t4ddg_1(0,2,0,0) - (t4ddg_2(0,2,0,0)*t2s_2(0,0))
		,"T2s(k,l)&T4ddg(i,j,k,l)(0,2,0,0)");
  test_for_zero(t4ddg_1(0,2,0,1) - (t4ddg_2(0,2,0,1)*t2s_2(0,1))
		,"T2s(k,l)&T4ddg(i,j,k,l)(0,2,0,1)");
  test_for_zero(t4ddg_1(0,2,0,2) - (t4ddg_2(0,2,0,2)*t2s_2(0,2))
		,"T2s(k,l)&T4ddg(i,j,k,l)(0,2,0,2)");
  test_for_zero(t4ddg_1(0,2,1,0) - (t4ddg_2(0,2,1,0)*t2s_2(1,0))
		,"T2s(k,l)&T4ddg(i,j,k,l)(0,2,1,0)");
  test_for_zero(t4ddg_1(0,2,1,1) - (t4ddg_2(0,2,1,1)*t2s_2(1,1))
		,"T2s(k,l)&T4ddg(i,j,k,l)(0,2,1,1)");
  test_for_zero(t4ddg_1(0,2,1,2) - (t4ddg_2(0,2,1,2)*t2s_2(1,2))
		,"T2s(k,l)&T4ddg(i,j,k,l)(0,2,1,2)");
  test_for_zero(t4ddg_1(0,2,2,0) - (t4ddg_2(0,2,2,0)*t2s_2(2,0))
		,"T2s(k,l)&T4ddg(i,j,k,l)(0,2,2,0)");
  test_for_zero(t4ddg_1(0,2,2,1) - (t4ddg_2(0,2,2,1)*t2s_2(2,1))
		,"T2s(k,l)&T4ddg(i,j,k,l)(0,2,2,1)");
  test_for_zero(t4ddg_1(0,2,2,2) - (t4ddg_2(0,2,2,2)*t2s_2(2,2))
		,"T2s(k,l)&T4ddg(i,j,k,l)(0,2,2,2)");
  test_for_zero(t4ddg_1(1,0,0,0) - (t4ddg_2(1,0,0,0)*t2s_2(0,0))
		,"T2s(k,l)&T4ddg(i,j,k,l)(1,0,0,0)");
  test_for_zero(t4ddg_1(1,0,0,1) - (t4ddg_2(1,0,0,1)*t2s_2(0,1))
		,"T2s(k,l)&T4ddg(i,j,k,l)(1,0,0,1)");
  test_for_zero(t4ddg_1(1,0,0,2) - (t4ddg_2(1,0,0,2)*t2s_2(0,2))
		,"T2s(k,l)&T4ddg(i,j,k,l)(1,0,0,2)");
  test_for_zero(t4ddg_1(1,0,1,0) - (t4ddg_2(1,0,1,0)*t2s_2(1,0))
		,"T2s(k,l)&T4ddg(i,j,k,l)(1,0,1,0)");
  test_for_zero(t4ddg_1(1,0,1,1) - (t4ddg_2(1,0,1,1)*t2s_2(1,1))
		,"T2s(k,l)&T4ddg(i,j,k,l)(1,0,1,1)");
  test_for_zero(t4ddg_1(1,0,1,2) - (t4ddg_2(1,0,1,2)*t2s_2(1,2))
		,"T2s(k,l)&T4ddg(i,j,k,l)(1,0,1,2)");
  test_for_zero(t4ddg_1(1,0,2,0) - (t4ddg_2(1,0,2,0)*t2s_2(2,0))
		,"T2s(k,l)&T4ddg(i,j,k,l)(1,0,2,0)");
  test_for_zero(t4ddg_1(1,0,2,1) - (t4ddg_2(1,0,2,1)*t2s_2(2,1))
		,"T2s(k,l)&T4ddg(i,j,k,l)(1,0,2,1)");
  test_for_zero(t4ddg_1(1,0,2,2) - (t4ddg_2(1,0,2,2)*t2s_2(2,2))
		,"T2s(k,l)&T4ddg(i,j,k,l)(1,0,2,2)");
  test_for_zero(t4ddg_1(1,1,0,0) - (t4ddg_2(1,1,0,0)*t2s_2(0,0))
		,"T2s(k,l)&T4ddg(i,j,k,l)(1,1,0,0)");
  test_for_zero(t4ddg_1(1,1,0,1) - (t4ddg_2(1,1,0,1)*t2s_2(0,1))
		,"T2s(k,l)&T4ddg(i,j,k,l)(1,1,0,1)");
  test_for_zero(t4ddg_1(1,1,0,2) - (t4ddg_2(1,1,0,2)*t2s_2(0,2))
		,"T2s(k,l)&T4ddg(i,j,k,l)(1,1,0,2)");
  test_for_zero(t4ddg_1(1,1,1,0) - (t4ddg_2(1,1,1,0)*t2s_2(1,0))
		,"T2s(k,l)&T4ddg(i,j,k,l)(1,1,1,0)");
  test_for_zero(t4ddg_1(1,1,1,1) - (t4ddg_2(1,1,1,1)*t2s_2(1,1))
		,"T2s(k,l)&T4ddg(i,j,k,l)(1,1,1,1)");
  test_for_zero(t4ddg_1(1,1,1,2) - (t4ddg_2(1,1,1,2)*t2s_2(1,2))
		,"T2s(k,l)&T4ddg(i,j,k,l)(1,1,1,2)");
  test_for_zero(t4ddg_1(1,1,2,0) - (t4ddg_2(1,1,2,0)*t2s_2(2,0))
		,"T2s(k,l)&T4ddg(i,j,k,l)(1,1,2,0)");
  test_for_zero(t4ddg_1(1,1,2,1) - (t4ddg_2(1,1,2,1)*t2s_2(2,1))
		,"T2s(k,l)&T4ddg(i,j,k,l)(1,1,2,1)");
  test_for_zero(t4ddg_1(1,1,2,2) - (t4ddg_2(1,1,2,2)*t2s_2(2,2))
		,"T2s(k,l)&T4ddg(i,j,k,l)(1,1,2,2)");
  test_for_zero(t4ddg_1(1,2,0,0) - (t4ddg_2(1,2,0,0)*t2s_2(0,0))
		,"T2s(k,l)&T4ddg(i,j,k,l)(1,2,0,0)");
  test_for_zero(t4ddg_1(1,2,0,1) - (t4ddg_2(1,2,0,1)*t2s_2(0,1))
		,"T2s(k,l)&T4ddg(i,j,k,l)(1,2,0,1)");
  test_for_zero(t4ddg_1(1,2,0,2) - (t4ddg_2(1,2,0,2)*t2s_2(0,2))
		,"T2s(k,l)&T4ddg(i,j,k,l)(1,2,0,2)");
  test_for_zero(t4ddg_1(1,2,1,0) - (t4ddg_2(1,2,1,0)*t2s_2(1,0))
		,"T2s(k,l)&T4ddg(i,j,k,l)(1,2,1,0)");
  test_for_zero(t4ddg_1(1,2,1,1) - (t4ddg_2(1,2,1,1)*t2s_2(1,1))
		,"T2s(k,l)&T4ddg(i,j,k,l)(1,2,1,1)");
  test_for_zero(t4ddg_1(1,2,1,2) - (t4ddg_2(1,2,1,2)*t2s_2(1,2))
		,"T2s(k,l)&T4ddg(i,j,k,l)(1,2,1,2)");
  test_for_zero(t4ddg_1(1,2,2,0) - (t4ddg_2(1,2,2,0)*t2s_2(2,0))
		,"T2s(k,l)&T4ddg(i,j,k,l)(1,2,2,0)");
  test_for_zero(t4ddg_1(1,2,2,1) - (t4ddg_2(1,2,2,1)*t2s_2(2,1))
		,"T2s(k,l)&T4ddg(i,j,k,l)(1,2,2,1)");
  test_for_zero(t4ddg_1(1,2,2,2) - (t4ddg_2(1,2,2,2)*t2s_2(2,2))
		,"T2s(k,l)&T4ddg(i,j,k,l)(1,2,2,2)");
  test_for_zero(t4ddg_1(2,0,0,0) - (t4ddg_2(2,0,0,0)*t2s_2(0,0))
		,"T2s(k,l)&T4ddg(i,j,k,l)(2,0,0,0)");
  test_for_zero(t4ddg_1(2,0,0,1) - (t4ddg_2(2,0,0,1)*t2s_2(0,1))
		,"T2s(k,l)&T4ddg(i,j,k,l)(2,0,0,1)");
  test_for_zero(t4ddg_1(2,0,0,2) - (t4ddg_2(2,0,0,2)*t2s_2(0,2))
		,"T2s(k,l)&T4ddg(i,j,k,l)(2,0,0,2)");
  test_for_zero(t4ddg_1(2,0,1,0) - (t4ddg_2(2,0,1,0)*t2s_2(1,0))
		,"T2s(k,l)&T4ddg(i,j,k,l)(2,0,1,0)");
  test_for_zero(t4ddg_1(2,0,1,1) - (t4ddg_2(2,0,1,1)*t2s_2(1,1))
		,"T2s(k,l)&T4ddg(i,j,k,l)(2,0,1,1)");
  test_for_zero(t4ddg_1(2,0,1,2) - (t4ddg_2(2,0,1,2)*t2s_2(1,2))
		,"T2s(k,l)&T4ddg(i,j,k,l)(2,0,1,2)");
  test_for_zero(t4ddg_1(2,0,2,0) - (t4ddg_2(2,0,2,0)*t2s_2(2,0))
		,"T2s(k,l)&T4ddg(i,j,k,l)(2,0,2,0)");
  test_for_zero(t4ddg_1(2,0,2,1) - (t4ddg_2(2,0,2,1)*t2s_2(2,1))
		,"T2s(k,l)&T4ddg(i,j,k,l)(2,0,2,1)");
  test_for_zero(t4ddg_1(2,0,2,2) - (t4ddg_2(2,0,2,2)*t2s_2(2,2))
		,"T2s(k,l)&T4ddg(i,j,k,l)(2,0,2,2)");
  test_for_zero(t4ddg_1(2,1,0,0) - (t4ddg_2(2,1,0,0)*t2s_2(0,0))
		,"T2s(k,l)&T4ddg(i,j,k,l)(2,1,0,0)");
  test_for_zero(t4ddg_1(2,1,0,1) - (t4ddg_2(2,1,0,1)*t2s_2(0,1))
		,"T2s(k,l)&T4ddg(i,j,k,l)(2,1,0,1)");
  test_for_zero(t4ddg_1(2,1,0,2) - (t4ddg_2(2,1,0,2)*t2s_2(0,2))
		,"T2s(k,l)&T4ddg(i,j,k,l)(2,1,0,2)");
  test_for_zero(t4ddg_1(2,1,1,0) - (t4ddg_2(2,1,1,0)*t2s_2(1,0))
		,"T2s(k,l)&T4ddg(i,j,k,l)(2,1,1,0)");
  test_for_zero(t4ddg_1(2,1,1,1) - (t4ddg_2(2,1,1,1)*t2s_2(1,1))
		,"T2s(k,l)&T4ddg(i,j,k,l)(2,1,1,1)");
  test_for_zero(t4ddg_1(2,1,1,2) - (t4ddg_2(2,1,1,2)*t2s_2(1,2))
		,"T2s(k,l)&T4ddg(i,j,k,l)(2,1,1,2)");
  test_for_zero(t4ddg_1(2,1,2,0) - (t4ddg_2(2,1,2,0)*t2s_2(2,0))
		,"T2s(k,l)&T4ddg(i,j,k,l)(2,1,2,0)");
  test_for_zero(t4ddg_1(2,1,2,1) - (t4ddg_2(2,1,2,1)*t2s_2(2,1))
		,"T2s(k,l)&T4ddg(i,j,k,l)(2,1,2,1)");
  test_for_zero(t4ddg_1(2,1,2,2) - (t4ddg_2(2,1,2,2)*t2s_2(2,2))
		,"T2s(k,l)&T4ddg(i,j,k,l)(2,1,2,2)");
  test_for_zero(t4ddg_1(2,2,0,0) - (t4ddg_2(2,2,0,0)*t2s_2(0,0))
		,"T2s(k,l)&T4ddg(i,j,k,l)(2,2,0,0)");
  test_for_zero(t4ddg_1(2,2,0,1) - (t4ddg_2(2,2,0,1)*t2s_2(0,1))
		,"T2s(k,l)&T4ddg(i,j,k,l)(2,2,0,1)");
  test_for_zero(t4ddg_1(2,2,0,2) - (t4ddg_2(2,2,0,2)*t2s_2(0,2))
		,"T2s(k,l)&T4ddg(i,j,k,l)(2,2,0,2)");
  test_for_zero(t4ddg_1(2,2,1,0) - (t4ddg_2(2,2,1,0)*t2s_2(1,0))
		,"T2s(k,l)&T4ddg(i,j,k,l)(2,2,1,0)");
  test_for_zero(t4ddg_1(2,2,1,1) - (t4ddg_2(2,2,1,1)*t2s_2(1,1))
		,"T2s(k,l)&T4ddg(i,j,k,l)(2,2,1,1)");
  test_for_zero(t4ddg_1(2,2,1,2) - (t4ddg_2(2,2,1,2)*t2s_2(1,2))
		,"T2s(k,l)&T4ddg(i,j,k,l)(2,2,1,2)");
  test_for_zero(t4ddg_1(2,2,2,0) - (t4ddg_2(2,2,2,0)*t2s_2(2,0))
		,"T2s(k,l)&T4ddg(i,j,k,l)(2,2,2,0)");
  test_for_zero(t4ddg_1(2,2,2,1) - (t4ddg_2(2,2,2,1)*t2s_2(2,1))
		,"T2s(k,l)&T4ddg(i,j,k,l)(2,2,2,1)");
  test_for_zero(t4ddg_1(2,2,2,2) - (t4ddg_2(2,2,2,2)*t2s_2(2,2))
		,"T2s(k,l)&T4ddg(i,j,k,l)(2,2,2,2)");

  t4ddg_1(i,j,k,l)=(t4ddg_2(i,j,k,l)%t2s_3(i,j));
  test_for_zero(t4ddg_1(0,0,0,0) - (t4ddg_2(0,0,0,0)/t2s_3(0,0))
		,"T4ddg(i,j,k,l)%T2s(i,j)(0,0,0,0)");
  test_for_zero(t4ddg_1(0,0,0,1) - (t4ddg_2(0,0,0,1)/t2s_3(0,0))
		,"T4ddg(i,j,k,l)%T2s(i,j)(0,0,0,1)");
  test_for_zero(t4ddg_1(0,0,0,2) - (t4ddg_2(0,0,0,2)/t2s_3(0,0))
		,"T4ddg(i,j,k,l)%T2s(i,j)(0,0,0,2)");
  test_for_zero(t4ddg_1(0,0,1,0) - (t4ddg_2(0,0,1,0)/t2s_3(0,0))
		,"T4ddg(i,j,k,l)%T2s(i,j)(0,0,1,0)");
  test_for_zero(t4ddg_1(0,0,1,1) - (t4ddg_2(0,0,1,1)/t2s_3(0,0))
		,"T4ddg(i,j,k,l)%T2s(i,j)(0,0,1,1)");
  test_for_zero(t4ddg_1(0,0,1,2) - (t4ddg_2(0,0,1,2)/t2s_3(0,0))
		,"T4ddg(i,j,k,l)%T2s(i,j)(0,0,1,2)");
  test_for_zero(t4ddg_1(0,0,2,0) - (t4ddg_2(0,0,2,0)/t2s_3(0,0))
		,"T4ddg(i,j,k,l)%T2s(i,j)(0,0,2,0)");
  test_for_zero(t4ddg_1(0,0,2,1) - (t4ddg_2(0,0,2,1)/t2s_3(0,0))
		,"T4ddg(i,j,k,l)%T2s(i,j)(0,0,2,1)");
  test_for_zero(t4ddg_1(0,0,2,2) - (t4ddg_2(0,0,2,2)/t2s_3(0,0))
		,"T4ddg(i,j,k,l)%T2s(i,j)(0,0,2,2)");
  test_for_zero(t4ddg_1(0,1,0,0) - (t4ddg_2(0,1,0,0)/t2s_3(0,1))
		,"T4ddg(i,j,k,l)%T2s(i,j)(0,1,0,0)");
  test_for_zero(t4ddg_1(0,1,0,1) - (t4ddg_2(0,1,0,1)/t2s_3(0,1))
		,"T4ddg(i,j,k,l)%T2s(i,j)(0,1,0,1)");
  test_for_zero(t4ddg_1(0,1,0,2) - (t4ddg_2(0,1,0,2)/t2s_3(0,1))
		,"T4ddg(i,j,k,l)%T2s(i,j)(0,1,0,2)");
  test_for_zero(t4ddg_1(0,1,1,0) - (t4ddg_2(0,1,1,0)/t2s_3(0,1))
		,"T4ddg(i,j,k,l)%T2s(i,j)(0,1,1,0)");
  test_for_zero(t4ddg_1(0,1,1,1) - (t4ddg_2(0,1,1,1)/t2s_3(0,1))
		,"T4ddg(i,j,k,l)%T2s(i,j)(0,1,1,1)");
  test_for_zero(t4ddg_1(0,1,1,2) - (t4ddg_2(0,1,1,2)/t2s_3(0,1))
		,"T4ddg(i,j,k,l)%T2s(i,j)(0,1,1,2)");
  test_for_zero(t4ddg_1(0,1,2,0) - (t4ddg_2(0,1,2,0)/t2s_3(0,1))
		,"T4ddg(i,j,k,l)%T2s(i,j)(0,1,2,0)");
  test_for_zero(t4ddg_1(0,1,2,1) - (t4ddg_2(0,1,2,1)/t2s_3(0,1))
		,"T4ddg(i,j,k,l)%T2s(i,j)(0,1,2,1)");
  test_for_zero(t4ddg_1(0,1,2,2) - (t4ddg_2(0,1,2,2)/t2s_3(0,1))
		,"T4ddg(i,j,k,l)%T2s(i,j)(0,1,2,2)");
  test_for_zero(t4ddg_1(0,2,0,0) - (t4ddg_2(0,2,0,0)/t2s_3(0,2))
		,"T4ddg(i,j,k,l)%T2s(i,j)(0,2,0,0)");
  test_for_zero(t4ddg_1(0,2,0,1) - (t4ddg_2(0,2,0,1)/t2s_3(0,2))
		,"T4ddg(i,j,k,l)%T2s(i,j)(0,2,0,1)");
  test_for_zero(t4ddg_1(0,2,0,2) - (t4ddg_2(0,2,0,2)/t2s_3(0,2))
		,"T4ddg(i,j,k,l)%T2s(i,j)(0,2,0,2)");
  test_for_zero(t4ddg_1(0,2,1,0) - (t4ddg_2(0,2,1,0)/t2s_3(0,2))
		,"T4ddg(i,j,k,l)%T2s(i,j)(0,2,1,0)");
  test_for_zero(t4ddg_1(0,2,1,1) - (t4ddg_2(0,2,1,1)/t2s_3(0,2))
		,"T4ddg(i,j,k,l)%T2s(i,j)(0,2,1,1)");
  test_for_zero(t4ddg_1(0,2,1,2) - (t4ddg_2(0,2,1,2)/t2s_3(0,2))
		,"T4ddg(i,j,k,l)%T2s(i,j)(0,2,1,2)");
  test_for_zero(t4ddg_1(0,2,2,0) - (t4ddg_2(0,2,2,0)/t2s_3(0,2))
		,"T4ddg(i,j,k,l)%T2s(i,j)(0,2,2,0)");
  test_for_zero(t4ddg_1(0,2,2,1) - (t4ddg_2(0,2,2,1)/t2s_3(0,2))
		,"T4ddg(i,j,k,l)%T2s(i,j)(0,2,2,1)");
  test_for_zero(t4ddg_1(0,2,2,2) - (t4ddg_2(0,2,2,2)/t2s_3(0,2))
		,"T4ddg(i,j,k,l)%T2s(i,j)(0,2,2,2)");
  test_for_zero(t4ddg_1(1,0,0,0) - (t4ddg_2(1,0,0,0)/t2s_3(1,0))
		,"T4ddg(i,j,k,l)%T2s(i,j)(1,0,0,0)");
  test_for_zero(t4ddg_1(1,0,0,1) - (t4ddg_2(1,0,0,1)/t2s_3(1,0))
		,"T4ddg(i,j,k,l)%T2s(i,j)(1,0,0,1)");
  test_for_zero(t4ddg_1(1,0,0,2) - (t4ddg_2(1,0,0,2)/t2s_3(1,0))
		,"T4ddg(i,j,k,l)%T2s(i,j)(1,0,0,2)");
  test_for_zero(t4ddg_1(1,0,1,0) - (t4ddg_2(1,0,1,0)/t2s_3(1,0))
		,"T4ddg(i,j,k,l)%T2s(i,j)(1,0,1,0)");
  test_for_zero(t4ddg_1(1,0,1,1) - (t4ddg_2(1,0,1,1)/t2s_3(1,0))
		,"T4ddg(i,j,k,l)%T2s(i,j)(1,0,1,1)");
  test_for_zero(t4ddg_1(1,0,1,2) - (t4ddg_2(1,0,1,2)/t2s_3(1,0))
		,"T4ddg(i,j,k,l)%T2s(i,j)(1,0,1,2)");
  test_for_zero(t4ddg_1(1,0,2,0) - (t4ddg_2(1,0,2,0)/t2s_3(1,0))
		,"T4ddg(i,j,k,l)%T2s(i,j)(1,0,2,0)");
  test_for_zero(t4ddg_1(1,0,2,1) - (t4ddg_2(1,0,2,1)/t2s_3(1,0))
		,"T4ddg(i,j,k,l)%T2s(i,j)(1,0,2,1)");
  test_for_zero(t4ddg_1(1,0,2,2) - (t4ddg_2(1,0,2,2)/t2s_3(1,0))
		,"T4ddg(i,j,k,l)%T2s(i,j)(1,0,2,2)");
  test_for_zero(t4ddg_1(1,1,0,0) - (t4ddg_2(1,1,0,0)/t2s_3(1,1))
		,"T4ddg(i,j,k,l)%T2s(i,j)(1,1,0,0)");
  test_for_zero(t4ddg_1(1,1,0,1) - (t4ddg_2(1,1,0,1)/t2s_3(1,1))
		,"T4ddg(i,j,k,l)%T2s(i,j)(1,1,0,1)");
  test_for_zero(t4ddg_1(1,1,0,2) - (t4ddg_2(1,1,0,2)/t2s_3(1,1))
		,"T4ddg(i,j,k,l)%T2s(i,j)(1,1,0,2)");
  test_for_zero(t4ddg_1(1,1,1,0) - (t4ddg_2(1,1,1,0)/t2s_3(1,1))
		,"T4ddg(i,j,k,l)%T2s(i,j)(1,1,1,0)");
  test_for_zero(t4ddg_1(1,1,1,1) - (t4ddg_2(1,1,1,1)/t2s_3(1,1))
		,"T4ddg(i,j,k,l)%T2s(i,j)(1,1,1,1)");
  test_for_zero(t4ddg_1(1,1,1,2) - (t4ddg_2(1,1,1,2)/t2s_3(1,1))
		,"T4ddg(i,j,k,l)%T2s(i,j)(1,1,1,2)");
  test_for_zero(t4ddg_1(1,1,2,0) - (t4ddg_2(1,1,2,0)/t2s_3(1,1))
		,"T4ddg(i,j,k,l)%T2s(i,j)(1,1,2,0)");
  test_for_zero(t4ddg_1(1,1,2,1) - (t4ddg_2(1,1,2,1)/t2s_3(1,1))
		,"T4ddg(i,j,k,l)%T2s(i,j)(1,1,2,1)");
  test_for_zero(t4ddg_1(1,1,2,2) - (t4ddg_2(1,1,2,2)/t2s_3(1,1))
		,"T4ddg(i,j,k,l)%T2s(i,j)(1,1,2,2)");
  test_for_zero(t4ddg_1(1,2,0,0) - (t4ddg_2(1,2,0,0)/t2s_3(1,2))
		,"T4ddg(i,j,k,l)%T2s(i,j)(1,2,0,0)");
  test_for_zero(t4ddg_1(1,2,0,1) - (t4ddg_2(1,2,0,1)/t2s_3(1,2))
		,"T4ddg(i,j,k,l)%T2s(i,j)(1,2,0,1)");
  test_for_zero(t4ddg_1(1,2,0,2) - (t4ddg_2(1,2,0,2)/t2s_3(1,2))
		,"T4ddg(i,j,k,l)%T2s(i,j)(1,2,0,2)");
  test_for_zero(t4ddg_1(1,2,1,0) - (t4ddg_2(1,2,1,0)/t2s_3(1,2))
		,"T4ddg(i,j,k,l)%T2s(i,j)(1,2,1,0)");
  test_for_zero(t4ddg_1(1,2,1,1) - (t4ddg_2(1,2,1,1)/t2s_3(1,2))
		,"T4ddg(i,j,k,l)%T2s(i,j)(1,2,1,1)");
  test_for_zero(t4ddg_1(1,2,1,2) - (t4ddg_2(1,2,1,2)/t2s_3(1,2))
		,"T4ddg(i,j,k,l)%T2s(i,j)(1,2,1,2)");
  test_for_zero(t4ddg_1(1,2,2,0) - (t4ddg_2(1,2,2,0)/t2s_3(1,2))
		,"T4ddg(i,j,k,l)%T2s(i,j)(1,2,2,0)");
  test_for_zero(t4ddg_1(1,2,2,1) - (t4ddg_2(1,2,2,1)/t2s_3(1,2))
		,"T4ddg(i,j,k,l)%T2s(i,j)(1,2,2,1)");
  test_for_zero(t4ddg_1(1,2,2,2) - (t4ddg_2(1,2,2,2)/t2s_3(1,2))
		,"T4ddg(i,j,k,l)%T2s(i,j)(1,2,2,2)");
  test_for_zero(t4ddg_1(2,0,0,0) - (t4ddg_2(2,0,0,0)/t2s_3(2,0))
		,"T4ddg(i,j,k,l)%T2s(i,j)(2,0,0,0)");
  test_for_zero(t4ddg_1(2,0,0,1) - (t4ddg_2(2,0,0,1)/t2s_3(2,0))
		,"T4ddg(i,j,k,l)%T2s(i,j)(2,0,0,1)");
  test_for_zero(t4ddg_1(2,0,0,2) - (t4ddg_2(2,0,0,2)/t2s_3(2,0))
		,"T4ddg(i,j,k,l)%T2s(i,j)(2,0,0,2)");
  test_for_zero(t4ddg_1(2,0,1,0) - (t4ddg_2(2,0,1,0)/t2s_3(2,0))
		,"T4ddg(i,j,k,l)%T2s(i,j)(2,0,1,0)");
  test_for_zero(t4ddg_1(2,0,1,1) - (t4ddg_2(2,0,1,1)/t2s_3(2,0))
		,"T4ddg(i,j,k,l)%T2s(i,j)(2,0,1,1)");
  test_for_zero(t4ddg_1(2,0,1,2) - (t4ddg_2(2,0,1,2)/t2s_3(2,0))
		,"T4ddg(i,j,k,l)%T2s(i,j)(2,0,1,2)");
  test_for_zero(t4ddg_1(2,0,2,0) - (t4ddg_2(2,0,2,0)/t2s_3(2,0))
		,"T4ddg(i,j,k,l)%T2s(i,j)(2,0,2,0)");
  test_for_zero(t4ddg_1(2,0,2,1) - (t4ddg_2(2,0,2,1)/t2s_3(2,0))
		,"T4ddg(i,j,k,l)%T2s(i,j)(2,0,2,1)");
  test_for_zero(t4ddg_1(2,0,2,2) - (t4ddg_2(2,0,2,2)/t2s_3(2,0))
		,"T4ddg(i,j,k,l)%T2s(i,j)(2,0,2,2)");
  test_for_zero(t4ddg_1(2,1,0,0) - (t4ddg_2(2,1,0,0)/t2s_3(2,1))
		,"T4ddg(i,j,k,l)%T2s(i,j)(2,1,0,0)");
  test_for_zero(t4ddg_1(2,1,0,1) - (t4ddg_2(2,1,0,1)/t2s_3(2,1))
		,"T4ddg(i,j,k,l)%T2s(i,j)(2,1,0,1)");
  test_for_zero(t4ddg_1(2,1,0,2) - (t4ddg_2(2,1,0,2)/t2s_3(2,1))
		,"T4ddg(i,j,k,l)%T2s(i,j)(2,1,0,2)");
  test_for_zero(t4ddg_1(2,1,1,0) - (t4ddg_2(2,1,1,0)/t2s_3(2,1))
		,"T4ddg(i,j,k,l)%T2s(i,j)(2,1,1,0)");
  test_for_zero(t4ddg_1(2,1,1,1) - (t4ddg_2(2,1,1,1)/t2s_3(2,1))
		,"T4ddg(i,j,k,l)%T2s(i,j)(2,1,1,1)");
  test_for_zero(t4ddg_1(2,1,1,2) - (t4ddg_2(2,1,1,2)/t2s_3(2,1))
		,"T4ddg(i,j,k,l)%T2s(i,j)(2,1,1,2)");
  test_for_zero(t4ddg_1(2,1,2,0) - (t4ddg_2(2,1,2,0)/t2s_3(2,1))
		,"T4ddg(i,j,k,l)%T2s(i,j)(2,1,2,0)");
  test_for_zero(t4ddg_1(2,1,2,1) - (t4ddg_2(2,1,2,1)/t2s_3(2,1))
		,"T4ddg(i,j,k,l)%T2s(i,j)(2,1,2,1)");
  test_for_zero(t4ddg_1(2,1,2,2) - (t4ddg_2(2,1,2,2)/t2s_3(2,1))
		,"T4ddg(i,j,k,l)%T2s(i,j)(2,1,2,2)");
  test_for_zero(t4ddg_1(2,2,0,0) - (t4ddg_2(2,2,0,0)/t2s_3(2,2))
		,"T4ddg(i,j,k,l)%T2s(i,j)(2,2,0,0)");
  test_for_zero(t4ddg_1(2,2,0,1) - (t4ddg_2(2,2,0,1)/t2s_3(2,2))
		,"T4ddg(i,j,k,l)%T2s(i,j)(2,2,0,1)");
  test_for_zero(t4ddg_1(2,2,0,2) - (t4ddg_2(2,2,0,2)/t2s_3(2,2))
		,"T4ddg(i,j,k,l)%T2s(i,j)(2,2,0,2)");
  test_for_zero(t4ddg_1(2,2,1,0) - (t4ddg_2(2,2,1,0)/t2s_3(2,2))
		,"T4ddg(i,j,k,l)%T2s(i,j)(2,2,1,0)");
  test_for_zero(t4ddg_1(2,2,1,1) - (t4ddg_2(2,2,1,1)/t2s_3(2,2))
		,"T4ddg(i,j,k,l)%T2s(i,j)(2,2,1,1)");
  test_for_zero(t4ddg_1(2,2,1,2) - (t4ddg_2(2,2,1,2)/t2s_3(2,2))
		,"T4ddg(i,j,k,l)%T2s(i,j)(2,2,1,2)");
  test_for_zero(t4ddg_1(2,2,2,0) - (t4ddg_2(2,2,2,0)/t2s_3(2,2))
		,"T4ddg(i,j,k,l)%T2s(i,j)(2,2,2,0)");
  test_for_zero(t4ddg_1(2,2,2,1) - (t4ddg_2(2,2,2,1)/t2s_3(2,2))
		,"T4ddg(i,j,k,l)%T2s(i,j)(2,2,2,1)");
  test_for_zero(t4ddg_1(2,2,2,2) - (t4ddg_2(2,2,2,2)/t2s_3(2,2))
		,"T4ddg(i,j,k,l)%T2s(i,j)(2,2,2,2)");

  t4ddg_1(i,j,k,l)=(t2s_2(i,j)%t4ddg_2(i,j,k,l));
  test_for_zero(t4ddg_1(0,0,0,0) - (t4ddg_2(0,0,0,0)/t2s_2(0,0))
		,"T2s(i,j)%T4ddg(i,j,k,l)(0,0,0,0)");
  test_for_zero(t4ddg_1(0,0,0,1) - (t4ddg_2(0,0,0,1)/t2s_2(0,0))
		,"T2s(i,j)%T4ddg(i,j,k,l)(0,0,0,1)");
  test_for_zero(t4ddg_1(0,0,0,2) - (t4ddg_2(0,0,0,2)/t2s_2(0,0))
		,"T2s(i,j)%T4ddg(i,j,k,l)(0,0,0,2)");
  test_for_zero(t4ddg_1(0,0,1,0) - (t4ddg_2(0,0,1,0)/t2s_2(0,0))
		,"T2s(i,j)%T4ddg(i,j,k,l)(0,0,1,0)");
  test_for_zero(t4ddg_1(0,0,1,1) - (t4ddg_2(0,0,1,1)/t2s_2(0,0))
		,"T2s(i,j)%T4ddg(i,j,k,l)(0,0,1,1)");
  test_for_zero(t4ddg_1(0,0,1,2) - (t4ddg_2(0,0,1,2)/t2s_2(0,0))
		,"T2s(i,j)%T4ddg(i,j,k,l)(0,0,1,2)");
  test_for_zero(t4ddg_1(0,0,2,0) - (t4ddg_2(0,0,2,0)/t2s_2(0,0))
		,"T2s(i,j)%T4ddg(i,j,k,l)(0,0,2,0)");
  test_for_zero(t4ddg_1(0,0,2,1) - (t4ddg_2(0,0,2,1)/t2s_2(0,0))
		,"T2s(i,j)%T4ddg(i,j,k,l)(0,0,2,1)");
  test_for_zero(t4ddg_1(0,0,2,2) - (t4ddg_2(0,0,2,2)/t2s_2(0,0))
		,"T2s(i,j)%T4ddg(i,j,k,l)(0,0,2,2)");
  test_for_zero(t4ddg_1(0,1,0,0) - (t4ddg_2(0,1,0,0)/t2s_2(0,1))
		,"T2s(i,j)%T4ddg(i,j,k,l)(0,1,0,0)");
  test_for_zero(t4ddg_1(0,1,0,1) - (t4ddg_2(0,1,0,1)/t2s_2(0,1))
		,"T2s(i,j)%T4ddg(i,j,k,l)(0,1,0,1)");
  test_for_zero(t4ddg_1(0,1,0,2) - (t4ddg_2(0,1,0,2)/t2s_2(0,1))
		,"T2s(i,j)%T4ddg(i,j,k,l)(0,1,0,2)");
  test_for_zero(t4ddg_1(0,1,1,0) - (t4ddg_2(0,1,1,0)/t2s_2(0,1))
		,"T2s(i,j)%T4ddg(i,j,k,l)(0,1,1,0)");
  test_for_zero(t4ddg_1(0,1,1,1) - (t4ddg_2(0,1,1,1)/t2s_2(0,1))
		,"T2s(i,j)%T4ddg(i,j,k,l)(0,1,1,1)");
  test_for_zero(t4ddg_1(0,1,1,2) - (t4ddg_2(0,1,1,2)/t2s_2(0,1))
		,"T2s(i,j)%T4ddg(i,j,k,l)(0,1,1,2)");
  test_for_zero(t4ddg_1(0,1,2,0) - (t4ddg_2(0,1,2,0)/t2s_2(0,1))
		,"T2s(i,j)%T4ddg(i,j,k,l)(0,1,2,0)");
  test_for_zero(t4ddg_1(0,1,2,1) - (t4ddg_2(0,1,2,1)/t2s_2(0,1))
		,"T2s(i,j)%T4ddg(i,j,k,l)(0,1,2,1)");
  test_for_zero(t4ddg_1(0,1,2,2) - (t4ddg_2(0,1,2,2)/t2s_2(0,1))
		,"T2s(i,j)%T4ddg(i,j,k,l)(0,1,2,2)");
  test_for_zero(t4ddg_1(0,2,0,0) - (t4ddg_2(0,2,0,0)/t2s_2(0,2))
		,"T2s(i,j)%T4ddg(i,j,k,l)(0,2,0,0)");
  test_for_zero(t4ddg_1(0,2,0,1) - (t4ddg_2(0,2,0,1)/t2s_2(0,2))
		,"T2s(i,j)%T4ddg(i,j,k,l)(0,2,0,1)");
  test_for_zero(t4ddg_1(0,2,0,2) - (t4ddg_2(0,2,0,2)/t2s_2(0,2))
		,"T2s(i,j)%T4ddg(i,j,k,l)(0,2,0,2)");
  test_for_zero(t4ddg_1(0,2,1,0) - (t4ddg_2(0,2,1,0)/t2s_2(0,2))
		,"T2s(i,j)%T4ddg(i,j,k,l)(0,2,1,0)");
  test_for_zero(t4ddg_1(0,2,1,1) - (t4ddg_2(0,2,1,1)/t2s_2(0,2))
		,"T2s(i,j)%T4ddg(i,j,k,l)(0,2,1,1)");
  test_for_zero(t4ddg_1(0,2,1,2) - (t4ddg_2(0,2,1,2)/t2s_2(0,2))
		,"T2s(i,j)%T4ddg(i,j,k,l)(0,2,1,2)");
  test_for_zero(t4ddg_1(0,2,2,0) - (t4ddg_2(0,2,2,0)/t2s_2(0,2))
		,"T2s(i,j)%T4ddg(i,j,k,l)(0,2,2,0)");
  test_for_zero(t4ddg_1(0,2,2,1) - (t4ddg_2(0,2,2,1)/t2s_2(0,2))
		,"T2s(i,j)%T4ddg(i,j,k,l)(0,2,2,1)");
  test_for_zero(t4ddg_1(0,2,2,2) - (t4ddg_2(0,2,2,2)/t2s_2(0,2))
		,"T2s(i,j)%T4ddg(i,j,k,l)(0,2,2,2)");
  test_for_zero(t4ddg_1(1,0,0,0) - (t4ddg_2(1,0,0,0)/t2s_2(1,0))
		,"T2s(i,j)%T4ddg(i,j,k,l)(1,0,0,0)");
  test_for_zero(t4ddg_1(1,0,0,1) - (t4ddg_2(1,0,0,1)/t2s_2(1,0))
		,"T2s(i,j)%T4ddg(i,j,k,l)(1,0,0,1)");
  test_for_zero(t4ddg_1(1,0,0,2) - (t4ddg_2(1,0,0,2)/t2s_2(1,0))
		,"T2s(i,j)%T4ddg(i,j,k,l)(1,0,0,2)");
  test_for_zero(t4ddg_1(1,0,1,0) - (t4ddg_2(1,0,1,0)/t2s_2(1,0))
		,"T2s(i,j)%T4ddg(i,j,k,l)(1,0,1,0)");
  test_for_zero(t4ddg_1(1,0,1,1) - (t4ddg_2(1,0,1,1)/t2s_2(1,0))
		,"T2s(i,j)%T4ddg(i,j,k,l)(1,0,1,1)");
  test_for_zero(t4ddg_1(1,0,1,2) - (t4ddg_2(1,0,1,2)/t2s_2(1,0))
		,"T2s(i,j)%T4ddg(i,j,k,l)(1,0,1,2)");
  test_for_zero(t4ddg_1(1,0,2,0) - (t4ddg_2(1,0,2,0)/t2s_2(1,0))
		,"T2s(i,j)%T4ddg(i,j,k,l)(1,0,2,0)");
  test_for_zero(t4ddg_1(1,0,2,1) - (t4ddg_2(1,0,2,1)/t2s_2(1,0))
		,"T2s(i,j)%T4ddg(i,j,k,l)(1,0,2,1)");
  test_for_zero(t4ddg_1(1,0,2,2) - (t4ddg_2(1,0,2,2)/t2s_2(1,0))
		,"T2s(i,j)%T4ddg(i,j,k,l)(1,0,2,2)");
  test_for_zero(t4ddg_1(1,1,0,0) - (t4ddg_2(1,1,0,0)/t2s_2(1,1))
		,"T2s(i,j)%T4ddg(i,j,k,l)(1,1,0,0)");
  test_for_zero(t4ddg_1(1,1,0,1) - (t4ddg_2(1,1,0,1)/t2s_2(1,1))
		,"T2s(i,j)%T4ddg(i,j,k,l)(1,1,0,1)");
  test_for_zero(t4ddg_1(1,1,0,2) - (t4ddg_2(1,1,0,2)/t2s_2(1,1))
		,"T2s(i,j)%T4ddg(i,j,k,l)(1,1,0,2)");
  test_for_zero(t4ddg_1(1,1,1,0) - (t4ddg_2(1,1,1,0)/t2s_2(1,1))
		,"T2s(i,j)%T4ddg(i,j,k,l)(1,1,1,0)");
  test_for_zero(t4ddg_1(1,1,1,1) - (t4ddg_2(1,1,1,1)/t2s_2(1,1))
		,"T2s(i,j)%T4ddg(i,j,k,l)(1,1,1,1)");
  test_for_zero(t4ddg_1(1,1,1,2) - (t4ddg_2(1,1,1,2)/t2s_2(1,1))
		,"T2s(i,j)%T4ddg(i,j,k,l)(1,1,1,2)");
  test_for_zero(t4ddg_1(1,1,2,0) - (t4ddg_2(1,1,2,0)/t2s_2(1,1))
		,"T2s(i,j)%T4ddg(i,j,k,l)(1,1,2,0)");
  test_for_zero(t4ddg_1(1,1,2,1) - (t4ddg_2(1,1,2,1)/t2s_2(1,1))
		,"T2s(i,j)%T4ddg(i,j,k,l)(1,1,2,1)");
  test_for_zero(t4ddg_1(1,1,2,2) - (t4ddg_2(1,1,2,2)/t2s_2(1,1))
		,"T2s(i,j)%T4ddg(i,j,k,l)(1,1,2,2)");
  test_for_zero(t4ddg_1(1,2,0,0) - (t4ddg_2(1,2,0,0)/t2s_2(1,2))
		,"T2s(i,j)%T4ddg(i,j,k,l)(1,2,0,0)");
  test_for_zero(t4ddg_1(1,2,0,1) - (t4ddg_2(1,2,0,1)/t2s_2(1,2))
		,"T2s(i,j)%T4ddg(i,j,k,l)(1,2,0,1)");
  test_for_zero(t4ddg_1(1,2,0,2) - (t4ddg_2(1,2,0,2)/t2s_2(1,2))
		,"T2s(i,j)%T4ddg(i,j,k,l)(1,2,0,2)");
  test_for_zero(t4ddg_1(1,2,1,0) - (t4ddg_2(1,2,1,0)/t2s_2(1,2))
		,"T2s(i,j)%T4ddg(i,j,k,l)(1,2,1,0)");
  test_for_zero(t4ddg_1(1,2,1,1) - (t4ddg_2(1,2,1,1)/t2s_2(1,2))
		,"T2s(i,j)%T4ddg(i,j,k,l)(1,2,1,1)");
  test_for_zero(t4ddg_1(1,2,1,2) - (t4ddg_2(1,2,1,2)/t2s_2(1,2))
		,"T2s(i,j)%T4ddg(i,j,k,l)(1,2,1,2)");
  test_for_zero(t4ddg_1(1,2,2,0) - (t4ddg_2(1,2,2,0)/t2s_2(1,2))
		,"T2s(i,j)%T4ddg(i,j,k,l)(1,2,2,0)");
  test_for_zero(t4ddg_1(1,2,2,1) - (t4ddg_2(1,2,2,1)/t2s_2(1,2))
		,"T2s(i,j)%T4ddg(i,j,k,l)(1,2,2,1)");
  test_for_zero(t4ddg_1(1,2,2,2) - (t4ddg_2(1,2,2,2)/t2s_2(1,2))
		,"T2s(i,j)%T4ddg(i,j,k,l)(1,2,2,2)");
  test_for_zero(t4ddg_1(2,0,0,0) - (t4ddg_2(2,0,0,0)/t2s_2(2,0))
		,"T2s(i,j)%T4ddg(i,j,k,l)(2,0,0,0)");
  test_for_zero(t4ddg_1(2,0,0,1) - (t4ddg_2(2,0,0,1)/t2s_2(2,0))
		,"T2s(i,j)%T4ddg(i,j,k,l)(2,0,0,1)");
  test_for_zero(t4ddg_1(2,0,0,2) - (t4ddg_2(2,0,0,2)/t2s_2(2,0))
		,"T2s(i,j)%T4ddg(i,j,k,l)(2,0,0,2)");
  test_for_zero(t4ddg_1(2,0,1,0) - (t4ddg_2(2,0,1,0)/t2s_2(2,0))
		,"T2s(i,j)%T4ddg(i,j,k,l)(2,0,1,0)");
  test_for_zero(t4ddg_1(2,0,1,1) - (t4ddg_2(2,0,1,1)/t2s_2(2,0))
		,"T2s(i,j)%T4ddg(i,j,k,l)(2,0,1,1)");
  test_for_zero(t4ddg_1(2,0,1,2) - (t4ddg_2(2,0,1,2)/t2s_2(2,0))
		,"T2s(i,j)%T4ddg(i,j,k,l)(2,0,1,2)");
  test_for_zero(t4ddg_1(2,0,2,0) - (t4ddg_2(2,0,2,0)/t2s_2(2,0))
		,"T2s(i,j)%T4ddg(i,j,k,l)(2,0,2,0)");
  test_for_zero(t4ddg_1(2,0,2,1) - (t4ddg_2(2,0,2,1)/t2s_2(2,0))
		,"T2s(i,j)%T4ddg(i,j,k,l)(2,0,2,1)");
  test_for_zero(t4ddg_1(2,0,2,2) - (t4ddg_2(2,0,2,2)/t2s_2(2,0))
		,"T2s(i,j)%T4ddg(i,j,k,l)(2,0,2,2)");
  test_for_zero(t4ddg_1(2,1,0,0) - (t4ddg_2(2,1,0,0)/t2s_2(2,1))
		,"T2s(i,j)%T4ddg(i,j,k,l)(2,1,0,0)");
  test_for_zero(t4ddg_1(2,1,0,1) - (t4ddg_2(2,1,0,1)/t2s_2(2,1))
		,"T2s(i,j)%T4ddg(i,j,k,l)(2,1,0,1)");
  test_for_zero(t4ddg_1(2,1,0,2) - (t4ddg_2(2,1,0,2)/t2s_2(2,1))
		,"T2s(i,j)%T4ddg(i,j,k,l)(2,1,0,2)");
  test_for_zero(t4ddg_1(2,1,1,0) - (t4ddg_2(2,1,1,0)/t2s_2(2,1))
		,"T2s(i,j)%T4ddg(i,j,k,l)(2,1,1,0)");
  test_for_zero(t4ddg_1(2,1,1,1) - (t4ddg_2(2,1,1,1)/t2s_2(2,1))
		,"T2s(i,j)%T4ddg(i,j,k,l)(2,1,1,1)");
  test_for_zero(t4ddg_1(2,1,1,2) - (t4ddg_2(2,1,1,2)/t2s_2(2,1))
		,"T2s(i,j)%T4ddg(i,j,k,l)(2,1,1,2)");
  test_for_zero(t4ddg_1(2,1,2,0) - (t4ddg_2(2,1,2,0)/t2s_2(2,1))
		,"T2s(i,j)%T4ddg(i,j,k,l)(2,1,2,0)");
  test_for_zero(t4ddg_1(2,1,2,1) - (t4ddg_2(2,1,2,1)/t2s_2(2,1))
		,"T2s(i,j)%T4ddg(i,j,k,l)(2,1,2,1)");
  test_for_zero(t4ddg_1(2,1,2,2) - (t4ddg_2(2,1,2,2)/t2s_2(2,1))
		,"T2s(i,j)%T4ddg(i,j,k,l)(2,1,2,2)");
  test_for_zero(t4ddg_1(2,2,0,0) - (t4ddg_2(2,2,0,0)/t2s_2(2,2))
		,"T2s(i,j)%T4ddg(i,j,k,l)(2,2,0,0)");
  test_for_zero(t4ddg_1(2,2,0,1) - (t4ddg_2(2,2,0,1)/t2s_2(2,2))
		,"T2s(i,j)%T4ddg(i,j,k,l)(2,2,0,1)");
  test_for_zero(t4ddg_1(2,2,0,2) - (t4ddg_2(2,2,0,2)/t2s_2(2,2))
		,"T2s(i,j)%T4ddg(i,j,k,l)(2,2,0,2)");
  test_for_zero(t4ddg_1(2,2,1,0) - (t4ddg_2(2,2,1,0)/t2s_2(2,2))
		,"T2s(i,j)%T4ddg(i,j,k,l)(2,2,1,0)");
  test_for_zero(t4ddg_1(2,2,1,1) - (t4ddg_2(2,2,1,1)/t2s_2(2,2))
		,"T2s(i,j)%T4ddg(i,j,k,l)(2,2,1,1)");
  test_for_zero(t4ddg_1(2,2,1,2) - (t4ddg_2(2,2,1,2)/t2s_2(2,2))
		,"T2s(i,j)%T4ddg(i,j,k,l)(2,2,1,2)");
  test_for_zero(t4ddg_1(2,2,2,0) - (t4ddg_2(2,2,2,0)/t2s_2(2,2))
		,"T2s(i,j)%T4ddg(i,j,k,l)(2,2,2,0)");
  test_for_zero(t4ddg_1(2,2,2,1) - (t4ddg_2(2,2,2,1)/t2s_2(2,2))
		,"T2s(i,j)%T4ddg(i,j,k,l)(2,2,2,1)");
  test_for_zero(t4ddg_1(2,2,2,2) - (t4ddg_2(2,2,2,2)/t2s_2(2,2))
		,"T2s(i,j)%T4ddg(i,j,k,l)(2,2,2,2)");

  t4ddg_1(i,j,k,l)=(t4ddg_2(i,j,k,l)%t2s_3(k,l));
  test_for_zero(t4ddg_1(0,0,0,0) - (t4ddg_2(0,0,0,0)/t2s_3(0,0))
		,"T4ddg(i,j,k,l)%T2s(k,l)(0,0,0,0)");
  test_for_zero(t4ddg_1(0,0,0,1) - (t4ddg_2(0,0,0,1)/t2s_3(0,1))
		,"T4ddg(i,j,k,l)%T2s(k,l)(0,0,0,1)");
  test_for_zero(t4ddg_1(0,0,0,2) - (t4ddg_2(0,0,0,2)/t2s_3(0,2))
		,"T4ddg(i,j,k,l)%T2s(k,l)(0,0,0,2)");
  test_for_zero(t4ddg_1(0,0,1,0) - (t4ddg_2(0,0,1,0)/t2s_3(1,0))
		,"T4ddg(i,j,k,l)%T2s(k,l)(0,0,1,0)");
  test_for_zero(t4ddg_1(0,0,1,1) - (t4ddg_2(0,0,1,1)/t2s_3(1,1))
		,"T4ddg(i,j,k,l)%T2s(k,l)(0,0,1,1)");
  test_for_zero(t4ddg_1(0,0,1,2) - (t4ddg_2(0,0,1,2)/t2s_3(1,2))
		,"T4ddg(i,j,k,l)%T2s(k,l)(0,0,1,2)");
  test_for_zero(t4ddg_1(0,0,2,0) - (t4ddg_2(0,0,2,0)/t2s_3(2,0))
		,"T4ddg(i,j,k,l)%T2s(k,l)(0,0,2,0)");
  test_for_zero(t4ddg_1(0,0,2,1) - (t4ddg_2(0,0,2,1)/t2s_3(2,1))
		,"T4ddg(i,j,k,l)%T2s(k,l)(0,0,2,1)");
  test_for_zero(t4ddg_1(0,0,2,2) - (t4ddg_2(0,0,2,2)/t2s_3(2,2))
		,"T4ddg(i,j,k,l)%T2s(k,l)(0,0,2,2)");
  test_for_zero(t4ddg_1(0,1,0,0) - (t4ddg_2(0,1,0,0)/t2s_3(0,0))
		,"T4ddg(i,j,k,l)%T2s(k,l)(0,1,0,0)");
  test_for_zero(t4ddg_1(0,1,0,1) - (t4ddg_2(0,1,0,1)/t2s_3(0,1))
		,"T4ddg(i,j,k,l)%T2s(k,l)(0,1,0,1)");
  test_for_zero(t4ddg_1(0,1,0,2) - (t4ddg_2(0,1,0,2)/t2s_3(0,2))
		,"T4ddg(i,j,k,l)%T2s(k,l)(0,1,0,2)");
  test_for_zero(t4ddg_1(0,1,1,0) - (t4ddg_2(0,1,1,0)/t2s_3(1,0))
		,"T4ddg(i,j,k,l)%T2s(k,l)(0,1,1,0)");
  test_for_zero(t4ddg_1(0,1,1,1) - (t4ddg_2(0,1,1,1)/t2s_3(1,1))
		,"T4ddg(i,j,k,l)%T2s(k,l)(0,1,1,1)");
  test_for_zero(t4ddg_1(0,1,1,2) - (t4ddg_2(0,1,1,2)/t2s_3(1,2))
		,"T4ddg(i,j,k,l)%T2s(k,l)(0,1,1,2)");
  test_for_zero(t4ddg_1(0,1,2,0) - (t4ddg_2(0,1,2,0)/t2s_3(2,0))
		,"T4ddg(i,j,k,l)%T2s(k,l)(0,1,2,0)");
  test_for_zero(t4ddg_1(0,1,2,1) - (t4ddg_2(0,1,2,1)/t2s_3(2,1))
		,"T4ddg(i,j,k,l)%T2s(k,l)(0,1,2,1)");
  test_for_zero(t4ddg_1(0,1,2,2) - (t4ddg_2(0,1,2,2)/t2s_3(2,2))
		,"T4ddg(i,j,k,l)%T2s(k,l)(0,1,2,2)");
  test_for_zero(t4ddg_1(0,2,0,0) - (t4ddg_2(0,2,0,0)/t2s_3(0,0))
		,"T4ddg(i,j,k,l)%T2s(k,l)(0,2,0,0)");
  test_for_zero(t4ddg_1(0,2,0,1) - (t4ddg_2(0,2,0,1)/t2s_3(0,1))
		,"T4ddg(i,j,k,l)%T2s(k,l)(0,2,0,1)");
  test_for_zero(t4ddg_1(0,2,0,2) - (t4ddg_2(0,2,0,2)/t2s_3(0,2))
		,"T4ddg(i,j,k,l)%T2s(k,l)(0,2,0,2)");
  test_for_zero(t4ddg_1(0,2,1,0) - (t4ddg_2(0,2,1,0)/t2s_3(1,0))
		,"T4ddg(i,j,k,l)%T2s(k,l)(0,2,1,0)");
  test_for_zero(t4ddg_1(0,2,1,1) - (t4ddg_2(0,2,1,1)/t2s_3(1,1))
		,"T4ddg(i,j,k,l)%T2s(k,l)(0,2,1,1)");
  test_for_zero(t4ddg_1(0,2,1,2) - (t4ddg_2(0,2,1,2)/t2s_3(1,2))
		,"T4ddg(i,j,k,l)%T2s(k,l)(0,2,1,2)");
  test_for_zero(t4ddg_1(0,2,2,0) - (t4ddg_2(0,2,2,0)/t2s_3(2,0))
		,"T4ddg(i,j,k,l)%T2s(k,l)(0,2,2,0)");
  test_for_zero(t4ddg_1(0,2,2,1) - (t4ddg_2(0,2,2,1)/t2s_3(2,1))
		,"T4ddg(i,j,k,l)%T2s(k,l)(0,2,2,1)");
  test_for_zero(t4ddg_1(0,2,2,2) - (t4ddg_2(0,2,2,2)/t2s_3(2,2))
		,"T4ddg(i,j,k,l)%T2s(k,l)(0,2,2,2)");
  test_for_zero(t4ddg_1(1,0,0,0) - (t4ddg_2(1,0,0,0)/t2s_3(0,0))
		,"T4ddg(i,j,k,l)%T2s(k,l)(1,0,0,0)");
  test_for_zero(t4ddg_1(1,0,0,1) - (t4ddg_2(1,0,0,1)/t2s_3(0,1))
		,"T4ddg(i,j,k,l)%T2s(k,l)(1,0,0,1)");
  test_for_zero(t4ddg_1(1,0,0,2) - (t4ddg_2(1,0,0,2)/t2s_3(0,2))
		,"T4ddg(i,j,k,l)%T2s(k,l)(1,0,0,2)");
  test_for_zero(t4ddg_1(1,0,1,0) - (t4ddg_2(1,0,1,0)/t2s_3(1,0))
		,"T4ddg(i,j,k,l)%T2s(k,l)(1,0,1,0)");
  test_for_zero(t4ddg_1(1,0,1,1) - (t4ddg_2(1,0,1,1)/t2s_3(1,1))
		,"T4ddg(i,j,k,l)%T2s(k,l)(1,0,1,1)");
  test_for_zero(t4ddg_1(1,0,1,2) - (t4ddg_2(1,0,1,2)/t2s_3(1,2))
		,"T4ddg(i,j,k,l)%T2s(k,l)(1,0,1,2)");
  test_for_zero(t4ddg_1(1,0,2,0) - (t4ddg_2(1,0,2,0)/t2s_3(2,0))
		,"T4ddg(i,j,k,l)%T2s(k,l)(1,0,2,0)");
  test_for_zero(t4ddg_1(1,0,2,1) - (t4ddg_2(1,0,2,1)/t2s_3(2,1))
		,"T4ddg(i,j,k,l)%T2s(k,l)(1,0,2,1)");
  test_for_zero(t4ddg_1(1,0,2,2) - (t4ddg_2(1,0,2,2)/t2s_3(2,2))
		,"T4ddg(i,j,k,l)%T2s(k,l)(1,0,2,2)");
  test_for_zero(t4ddg_1(1,1,0,0) - (t4ddg_2(1,1,0,0)/t2s_3(0,0))
		,"T4ddg(i,j,k,l)%T2s(k,l)(1,1,0,0)");
  test_for_zero(t4ddg_1(1,1,0,1) - (t4ddg_2(1,1,0,1)/t2s_3(0,1))
		,"T4ddg(i,j,k,l)%T2s(k,l)(1,1,0,1)");
  test_for_zero(t4ddg_1(1,1,0,2) - (t4ddg_2(1,1,0,2)/t2s_3(0,2))
		,"T4ddg(i,j,k,l)%T2s(k,l)(1,1,0,2)");
  test_for_zero(t4ddg_1(1,1,1,0) - (t4ddg_2(1,1,1,0)/t2s_3(1,0))
		,"T4ddg(i,j,k,l)%T2s(k,l)(1,1,1,0)");
  test_for_zero(t4ddg_1(1,1,1,1) - (t4ddg_2(1,1,1,1)/t2s_3(1,1))
		,"T4ddg(i,j,k,l)%T2s(k,l)(1,1,1,1)");
  test_for_zero(t4ddg_1(1,1,1,2) - (t4ddg_2(1,1,1,2)/t2s_3(1,2))
		,"T4ddg(i,j,k,l)%T2s(k,l)(1,1,1,2)");
  test_for_zero(t4ddg_1(1,1,2,0) - (t4ddg_2(1,1,2,0)/t2s_3(2,0))
		,"T4ddg(i,j,k,l)%T2s(k,l)(1,1,2,0)");
  test_for_zero(t4ddg_1(1,1,2,1) - (t4ddg_2(1,1,2,1)/t2s_3(2,1))
		,"T4ddg(i,j,k,l)%T2s(k,l)(1,1,2,1)");
  test_for_zero(t4ddg_1(1,1,2,2) - (t4ddg_2(1,1,2,2)/t2s_3(2,2))
		,"T4ddg(i,j,k,l)%T2s(k,l)(1,1,2,2)");
  test_for_zero(t4ddg_1(1,2,0,0) - (t4ddg_2(1,2,0,0)/t2s_3(0,0))
		,"T4ddg(i,j,k,l)%T2s(k,l)(1,2,0,0)");
  test_for_zero(t4ddg_1(1,2,0,1) - (t4ddg_2(1,2,0,1)/t2s_3(0,1))
		,"T4ddg(i,j,k,l)%T2s(k,l)(1,2,0,1)");
  test_for_zero(t4ddg_1(1,2,0,2) - (t4ddg_2(1,2,0,2)/t2s_3(0,2))
		,"T4ddg(i,j,k,l)%T2s(k,l)(1,2,0,2)");
  test_for_zero(t4ddg_1(1,2,1,0) - (t4ddg_2(1,2,1,0)/t2s_3(1,0))
		,"T4ddg(i,j,k,l)%T2s(k,l)(1,2,1,0)");
  test_for_zero(t4ddg_1(1,2,1,1) - (t4ddg_2(1,2,1,1)/t2s_3(1,1))
		,"T4ddg(i,j,k,l)%T2s(k,l)(1,2,1,1)");
  test_for_zero(t4ddg_1(1,2,1,2) - (t4ddg_2(1,2,1,2)/t2s_3(1,2))
		,"T4ddg(i,j,k,l)%T2s(k,l)(1,2,1,2)");
  test_for_zero(t4ddg_1(1,2,2,0) - (t4ddg_2(1,2,2,0)/t2s_3(2,0))
		,"T4ddg(i,j,k,l)%T2s(k,l)(1,2,2,0)");
  test_for_zero(t4ddg_1(1,2,2,1) - (t4ddg_2(1,2,2,1)/t2s_3(2,1))
		,"T4ddg(i,j,k,l)%T2s(k,l)(1,2,2,1)");
  test_for_zero(t4ddg_1(1,2,2,2) - (t4ddg_2(1,2,2,2)/t2s_3(2,2))
		,"T4ddg(i,j,k,l)%T2s(k,l)(1,2,2,2)");
  test_for_zero(t4ddg_1(2,0,0,0) - (t4ddg_2(2,0,0,0)/t2s_3(0,0))
		,"T4ddg(i,j,k,l)%T2s(k,l)(2,0,0,0)");
  test_for_zero(t4ddg_1(2,0,0,1) - (t4ddg_2(2,0,0,1)/t2s_3(0,1))
		,"T4ddg(i,j,k,l)%T2s(k,l)(2,0,0,1)");
  test_for_zero(t4ddg_1(2,0,0,2) - (t4ddg_2(2,0,0,2)/t2s_3(0,2))
		,"T4ddg(i,j,k,l)%T2s(k,l)(2,0,0,2)");
  test_for_zero(t4ddg_1(2,0,1,0) - (t4ddg_2(2,0,1,0)/t2s_3(1,0))
		,"T4ddg(i,j,k,l)%T2s(k,l)(2,0,1,0)");
  test_for_zero(t4ddg_1(2,0,1,1) - (t4ddg_2(2,0,1,1)/t2s_3(1,1))
		,"T4ddg(i,j,k,l)%T2s(k,l)(2,0,1,1)");
  test_for_zero(t4ddg_1(2,0,1,2) - (t4ddg_2(2,0,1,2)/t2s_3(1,2))
		,"T4ddg(i,j,k,l)%T2s(k,l)(2,0,1,2)");
  test_for_zero(t4ddg_1(2,0,2,0) - (t4ddg_2(2,0,2,0)/t2s_3(2,0))
		,"T4ddg(i,j,k,l)%T2s(k,l)(2,0,2,0)");
  test_for_zero(t4ddg_1(2,0,2,1) - (t4ddg_2(2,0,2,1)/t2s_3(2,1))
		,"T4ddg(i,j,k,l)%T2s(k,l)(2,0,2,1)");
  test_for_zero(t4ddg_1(2,0,2,2) - (t4ddg_2(2,0,2,2)/t2s_3(2,2))
		,"T4ddg(i,j,k,l)%T2s(k,l)(2,0,2,2)");
  test_for_zero(t4ddg_1(2,1,0,0) - (t4ddg_2(2,1,0,0)/t2s_3(0,0))
		,"T4ddg(i,j,k,l)%T2s(k,l)(2,1,0,0)");
  test_for_zero(t4ddg_1(2,1,0,1) - (t4ddg_2(2,1,0,1)/t2s_3(0,1))
		,"T4ddg(i,j,k,l)%T2s(k,l)(2,1,0,1)");
  test_for_zero(t4ddg_1(2,1,0,2) - (t4ddg_2(2,1,0,2)/t2s_3(0,2))
		,"T4ddg(i,j,k,l)%T2s(k,l)(2,1,0,2)");
  test_for_zero(t4ddg_1(2,1,1,0) - (t4ddg_2(2,1,1,0)/t2s_3(1,0))
		,"T4ddg(i,j,k,l)%T2s(k,l)(2,1,1,0)");
  test_for_zero(t4ddg_1(2,1,1,1) - (t4ddg_2(2,1,1,1)/t2s_3(1,1))
		,"T4ddg(i,j,k,l)%T2s(k,l)(2,1,1,1)");
  test_for_zero(t4ddg_1(2,1,1,2) - (t4ddg_2(2,1,1,2)/t2s_3(1,2))
		,"T4ddg(i,j,k,l)%T2s(k,l)(2,1,1,2)");
  test_for_zero(t4ddg_1(2,1,2,0) - (t4ddg_2(2,1,2,0)/t2s_3(2,0))
		,"T4ddg(i,j,k,l)%T2s(k,l)(2,1,2,0)");
  test_for_zero(t4ddg_1(2,1,2,1) - (t4ddg_2(2,1,2,1)/t2s_3(2,1))
		,"T4ddg(i,j,k,l)%T2s(k,l)(2,1,2,1)");
  test_for_zero(t4ddg_1(2,1,2,2) - (t4ddg_2(2,1,2,2)/t2s_3(2,2))
		,"T4ddg(i,j,k,l)%T2s(k,l)(2,1,2,2)");
  test_for_zero(t4ddg_1(2,2,0,0) - (t4ddg_2(2,2,0,0)/t2s_3(0,0))
		,"T4ddg(i,j,k,l)%T2s(k,l)(2,2,0,0)");
  test_for_zero(t4ddg_1(2,2,0,1) - (t4ddg_2(2,2,0,1)/t2s_3(0,1))
		,"T4ddg(i,j,k,l)%T2s(k,l)(2,2,0,1)");
  test_for_zero(t4ddg_1(2,2,0,2) - (t4ddg_2(2,2,0,2)/t2s_3(0,2))
		,"T4ddg(i,j,k,l)%T2s(k,l)(2,2,0,2)");
  test_for_zero(t4ddg_1(2,2,1,0) - (t4ddg_2(2,2,1,0)/t2s_3(1,0))
		,"T4ddg(i,j,k,l)%T2s(k,l)(2,2,1,0)");
  test_for_zero(t4ddg_1(2,2,1,1) - (t4ddg_2(2,2,1,1)/t2s_3(1,1))
		,"T4ddg(i,j,k,l)%T2s(k,l)(2,2,1,1)");
  test_for_zero(t4ddg_1(2,2,1,2) - (t4ddg_2(2,2,1,2)/t2s_3(1,2))
		,"T4ddg(i,j,k,l)%T2s(k,l)(2,2,1,2)");
  test_for_zero(t4ddg_1(2,2,2,0) - (t4ddg_2(2,2,2,0)/t2s_3(2,0))
		,"T4ddg(i,j,k,l)%T2s(k,l)(2,2,2,0)");
  test_for_zero(t4ddg_1(2,2,2,1) - (t4ddg_2(2,2,2,1)/t2s_3(2,1))
		,"T4ddg(i,j,k,l)%T2s(k,l)(2,2,2,1)");
  test_for_zero(t4ddg_1(2,2,2,2) - (t4ddg_2(2,2,2,2)/t2s_3(2,2))
		,"T4ddg(i,j,k,l)%T2s(k,l)(2,2,2,2)");

  t4ddg_1(i,j,k,l)=(t2s_2(k,l)%t4ddg_2(i,j,k,l));
  test_for_zero(t4ddg_1(0,0,0,0) - (t4ddg_2(0,0,0,0)/t2s_2(0,0))
		,"T2s(k,l)%T4ddg(i,j,k,l)(0,0,0,0)");
  test_for_zero(t4ddg_1(0,0,0,1) - (t4ddg_2(0,0,0,1)/t2s_2(0,1))
		,"T2s(k,l)%T4ddg(i,j,k,l)(0,0,0,1)");
  test_for_zero(t4ddg_1(0,0,0,2) - (t4ddg_2(0,0,0,2)/t2s_2(0,2))
		,"T2s(k,l)%T4ddg(i,j,k,l)(0,0,0,2)");
  test_for_zero(t4ddg_1(0,0,1,0) - (t4ddg_2(0,0,1,0)/t2s_2(1,0))
		,"T2s(k,l)%T4ddg(i,j,k,l)(0,0,1,0)");
  test_for_zero(t4ddg_1(0,0,1,1) - (t4ddg_2(0,0,1,1)/t2s_2(1,1))
		,"T2s(k,l)%T4ddg(i,j,k,l)(0,0,1,1)");
  test_for_zero(t4ddg_1(0,0,1,2) - (t4ddg_2(0,0,1,2)/t2s_2(1,2))
		,"T2s(k,l)%T4ddg(i,j,k,l)(0,0,1,2)");
  test_for_zero(t4ddg_1(0,0,2,0) - (t4ddg_2(0,0,2,0)/t2s_2(2,0))
		,"T2s(k,l)%T4ddg(i,j,k,l)(0,0,2,0)");
  test_for_zero(t4ddg_1(0,0,2,1) - (t4ddg_2(0,0,2,1)/t2s_2(2,1))
		,"T2s(k,l)%T4ddg(i,j,k,l)(0,0,2,1)");
  test_for_zero(t4ddg_1(0,0,2,2) - (t4ddg_2(0,0,2,2)/t2s_2(2,2))
		,"T2s(k,l)%T4ddg(i,j,k,l)(0,0,2,2)");
  test_for_zero(t4ddg_1(0,1,0,0) - (t4ddg_2(0,1,0,0)/t2s_2(0,0))
		,"T2s(k,l)%T4ddg(i,j,k,l)(0,1,0,0)");
  test_for_zero(t4ddg_1(0,1,0,1) - (t4ddg_2(0,1,0,1)/t2s_2(0,1))
		,"T2s(k,l)%T4ddg(i,j,k,l)(0,1,0,1)");
  test_for_zero(t4ddg_1(0,1,0,2) - (t4ddg_2(0,1,0,2)/t2s_2(0,2))
		,"T2s(k,l)%T4ddg(i,j,k,l)(0,1,0,2)");
  test_for_zero(t4ddg_1(0,1,1,0) - (t4ddg_2(0,1,1,0)/t2s_2(1,0))
		,"T2s(k,l)%T4ddg(i,j,k,l)(0,1,1,0)");
  test_for_zero(t4ddg_1(0,1,1,1) - (t4ddg_2(0,1,1,1)/t2s_2(1,1))
		,"T2s(k,l)%T4ddg(i,j,k,l)(0,1,1,1)");
  test_for_zero(t4ddg_1(0,1,1,2) - (t4ddg_2(0,1,1,2)/t2s_2(1,2))
		,"T2s(k,l)%T4ddg(i,j,k,l)(0,1,1,2)");
  test_for_zero(t4ddg_1(0,1,2,0) - (t4ddg_2(0,1,2,0)/t2s_2(2,0))
		,"T2s(k,l)%T4ddg(i,j,k,l)(0,1,2,0)");
  test_for_zero(t4ddg_1(0,1,2,1) - (t4ddg_2(0,1,2,1)/t2s_2(2,1))
		,"T2s(k,l)%T4ddg(i,j,k,l)(0,1,2,1)");
  test_for_zero(t4ddg_1(0,1,2,2) - (t4ddg_2(0,1,2,2)/t2s_2(2,2))
		,"T2s(k,l)%T4ddg(i,j,k,l)(0,1,2,2)");
  test_for_zero(t4ddg_1(0,2,0,0) - (t4ddg_2(0,2,0,0)/t2s_2(0,0))
		,"T2s(k,l)%T4ddg(i,j,k,l)(0,2,0,0)");
  test_for_zero(t4ddg_1(0,2,0,1) - (t4ddg_2(0,2,0,1)/t2s_2(0,1))
		,"T2s(k,l)%T4ddg(i,j,k,l)(0,2,0,1)");
  test_for_zero(t4ddg_1(0,2,0,2) - (t4ddg_2(0,2,0,2)/t2s_2(0,2))
		,"T2s(k,l)%T4ddg(i,j,k,l)(0,2,0,2)");
  test_for_zero(t4ddg_1(0,2,1,0) - (t4ddg_2(0,2,1,0)/t2s_2(1,0))
		,"T2s(k,l)%T4ddg(i,j,k,l)(0,2,1,0)");
  test_for_zero(t4ddg_1(0,2,1,1) - (t4ddg_2(0,2,1,1)/t2s_2(1,1))
		,"T2s(k,l)%T4ddg(i,j,k,l)(0,2,1,1)");
  test_for_zero(t4ddg_1(0,2,1,2) - (t4ddg_2(0,2,1,2)/t2s_2(1,2))
		,"T2s(k,l)%T4ddg(i,j,k,l)(0,2,1,2)");
  test_for_zero(t4ddg_1(0,2,2,0) - (t4ddg_2(0,2,2,0)/t2s_2(2,0))
		,"T2s(k,l)%T4ddg(i,j,k,l)(0,2,2,0)");
  test_for_zero(t4ddg_1(0,2,2,1) - (t4ddg_2(0,2,2,1)/t2s_2(2,1))
		,"T2s(k,l)%T4ddg(i,j,k,l)(0,2,2,1)");
  test_for_zero(t4ddg_1(0,2,2,2) - (t4ddg_2(0,2,2,2)/t2s_2(2,2))
		,"T2s(k,l)%T4ddg(i,j,k,l)(0,2,2,2)");
  test_for_zero(t4ddg_1(1,0,0,0) - (t4ddg_2(1,0,0,0)/t2s_2(0,0))
		,"T2s(k,l)%T4ddg(i,j,k,l)(1,0,0,0)");
  test_for_zero(t4ddg_1(1,0,0,1) - (t4ddg_2(1,0,0,1)/t2s_2(0,1))
		,"T2s(k,l)%T4ddg(i,j,k,l)(1,0,0,1)");
  test_for_zero(t4ddg_1(1,0,0,2) - (t4ddg_2(1,0,0,2)/t2s_2(0,2))
		,"T2s(k,l)%T4ddg(i,j,k,l)(1,0,0,2)");
  test_for_zero(t4ddg_1(1,0,1,0) - (t4ddg_2(1,0,1,0)/t2s_2(1,0))
		,"T2s(k,l)%T4ddg(i,j,k,l)(1,0,1,0)");
  test_for_zero(t4ddg_1(1,0,1,1) - (t4ddg_2(1,0,1,1)/t2s_2(1,1))
		,"T2s(k,l)%T4ddg(i,j,k,l)(1,0,1,1)");
  test_for_zero(t4ddg_1(1,0,1,2) - (t4ddg_2(1,0,1,2)/t2s_2(1,2))
		,"T2s(k,l)%T4ddg(i,j,k,l)(1,0,1,2)");
  test_for_zero(t4ddg_1(1,0,2,0) - (t4ddg_2(1,0,2,0)/t2s_2(2,0))
		,"T2s(k,l)%T4ddg(i,j,k,l)(1,0,2,0)");
  test_for_zero(t4ddg_1(1,0,2,1) - (t4ddg_2(1,0,2,1)/t2s_2(2,1))
		,"T2s(k,l)%T4ddg(i,j,k,l)(1,0,2,1)");
  test_for_zero(t4ddg_1(1,0,2,2) - (t4ddg_2(1,0,2,2)/t2s_2(2,2))
		,"T2s(k,l)%T4ddg(i,j,k,l)(1,0,2,2)");
  test_for_zero(t4ddg_1(1,1,0,0) - (t4ddg_2(1,1,0,0)/t2s_2(0,0))
		,"T2s(k,l)%T4ddg(i,j,k,l)(1,1,0,0)");
  test_for_zero(t4ddg_1(1,1,0,1) - (t4ddg_2(1,1,0,1)/t2s_2(0,1))
		,"T2s(k,l)%T4ddg(i,j,k,l)(1,1,0,1)");
  test_for_zero(t4ddg_1(1,1,0,2) - (t4ddg_2(1,1,0,2)/t2s_2(0,2))
		,"T2s(k,l)%T4ddg(i,j,k,l)(1,1,0,2)");
  test_for_zero(t4ddg_1(1,1,1,0) - (t4ddg_2(1,1,1,0)/t2s_2(1,0))
		,"T2s(k,l)%T4ddg(i,j,k,l)(1,1,1,0)");
  test_for_zero(t4ddg_1(1,1,1,1) - (t4ddg_2(1,1,1,1)/t2s_2(1,1))
		,"T2s(k,l)%T4ddg(i,j,k,l)(1,1,1,1)");
  test_for_zero(t4ddg_1(1,1,1,2) - (t4ddg_2(1,1,1,2)/t2s_2(1,2))
		,"T2s(k,l)%T4ddg(i,j,k,l)(1,1,1,2)");
  test_for_zero(t4ddg_1(1,1,2,0) - (t4ddg_2(1,1,2,0)/t2s_2(2,0))
		,"T2s(k,l)%T4ddg(i,j,k,l)(1,1,2,0)");
  test_for_zero(t4ddg_1(1,1,2,1) - (t4ddg_2(1,1,2,1)/t2s_2(2,1))
		,"T2s(k,l)%T4ddg(i,j,k,l)(1,1,2,1)");
  test_for_zero(t4ddg_1(1,1,2,2) - (t4ddg_2(1,1,2,2)/t2s_2(2,2))
		,"T2s(k,l)%T4ddg(i,j,k,l)(1,1,2,2)");
  test_for_zero(t4ddg_1(1,2,0,0) - (t4ddg_2(1,2,0,0)/t2s_2(0,0))
		,"T2s(k,l)%T4ddg(i,j,k,l)(1,2,0,0)");
  test_for_zero(t4ddg_1(1,2,0,1) - (t4ddg_2(1,2,0,1)/t2s_2(0,1))
		,"T2s(k,l)%T4ddg(i,j,k,l)(1,2,0,1)");
  test_for_zero(t4ddg_1(1,2,0,2) - (t4ddg_2(1,2,0,2)/t2s_2(0,2))
		,"T2s(k,l)%T4ddg(i,j,k,l)(1,2,0,2)");
  test_for_zero(t4ddg_1(1,2,1,0) - (t4ddg_2(1,2,1,0)/t2s_2(1,0))
		,"T2s(k,l)%T4ddg(i,j,k,l)(1,2,1,0)");
  test_for_zero(t4ddg_1(1,2,1,1) - (t4ddg_2(1,2,1,1)/t2s_2(1,1))
		,"T2s(k,l)%T4ddg(i,j,k,l)(1,2,1,1)");
  test_for_zero(t4ddg_1(1,2,1,2) - (t4ddg_2(1,2,1,2)/t2s_2(1,2))
		,"T2s(k,l)%T4ddg(i,j,k,l)(1,2,1,2)");
  test_for_zero(t4ddg_1(1,2,2,0) - (t4ddg_2(1,2,2,0)/t2s_2(2,0))
		,"T2s(k,l)%T4ddg(i,j,k,l)(1,2,2,0)");
  test_for_zero(t4ddg_1(1,2,2,1) - (t4ddg_2(1,2,2,1)/t2s_2(2,1))
		,"T2s(k,l)%T4ddg(i,j,k,l)(1,2,2,1)");
  test_for_zero(t4ddg_1(1,2,2,2) - (t4ddg_2(1,2,2,2)/t2s_2(2,2))
		,"T2s(k,l)%T4ddg(i,j,k,l)(1,2,2,2)");
  test_for_zero(t4ddg_1(2,0,0,0) - (t4ddg_2(2,0,0,0)/t2s_2(0,0))
		,"T2s(k,l)%T4ddg(i,j,k,l)(2,0,0,0)");
  test_for_zero(t4ddg_1(2,0,0,1) - (t4ddg_2(2,0,0,1)/t2s_2(0,1))
		,"T2s(k,l)%T4ddg(i,j,k,l)(2,0,0,1)");
  test_for_zero(t4ddg_1(2,0,0,2) - (t4ddg_2(2,0,0,2)/t2s_2(0,2))
		,"T2s(k,l)%T4ddg(i,j,k,l)(2,0,0,2)");
  test_for_zero(t4ddg_1(2,0,1,0) - (t4ddg_2(2,0,1,0)/t2s_2(1,0))
		,"T2s(k,l)%T4ddg(i,j,k,l)(2,0,1,0)");
  test_for_zero(t4ddg_1(2,0,1,1) - (t4ddg_2(2,0,1,1)/t2s_2(1,1))
		,"T2s(k,l)%T4ddg(i,j,k,l)(2,0,1,1)");
  test_for_zero(t4ddg_1(2,0,1,2) - (t4ddg_2(2,0,1,2)/t2s_2(1,2))
		,"T2s(k,l)%T4ddg(i,j,k,l)(2,0,1,2)");
  test_for_zero(t4ddg_1(2,0,2,0) - (t4ddg_2(2,0,2,0)/t2s_2(2,0))
		,"T2s(k,l)%T4ddg(i,j,k,l)(2,0,2,0)");
  test_for_zero(t4ddg_1(2,0,2,1) - (t4ddg_2(2,0,2,1)/t2s_2(2,1))
		,"T2s(k,l)%T4ddg(i,j,k,l)(2,0,2,1)");
  test_for_zero(t4ddg_1(2,0,2,2) - (t4ddg_2(2,0,2,2)/t2s_2(2,2))
		,"T2s(k,l)%T4ddg(i,j,k,l)(2,0,2,2)");
  test_for_zero(t4ddg_1(2,1,0,0) - (t4ddg_2(2,1,0,0)/t2s_2(0,0))
		,"T2s(k,l)%T4ddg(i,j,k,l)(2,1,0,0)");
  test_for_zero(t4ddg_1(2,1,0,1) - (t4ddg_2(2,1,0,1)/t2s_2(0,1))
		,"T2s(k,l)%T4ddg(i,j,k,l)(2,1,0,1)");
  test_for_zero(t4ddg_1(2,1,0,2) - (t4ddg_2(2,1,0,2)/t2s_2(0,2))
		,"T2s(k,l)%T4ddg(i,j,k,l)(2,1,0,2)");
  test_for_zero(t4ddg_1(2,1,1,0) - (t4ddg_2(2,1,1,0)/t2s_2(1,0))
		,"T2s(k,l)%T4ddg(i,j,k,l)(2,1,1,0)");
  test_for_zero(t4ddg_1(2,1,1,1) - (t4ddg_2(2,1,1,1)/t2s_2(1,1))
		,"T2s(k,l)%T4ddg(i,j,k,l)(2,1,1,1)");
  test_for_zero(t4ddg_1(2,1,1,2) - (t4ddg_2(2,1,1,2)/t2s_2(1,2))
		,"T2s(k,l)%T4ddg(i,j,k,l)(2,1,1,2)");
  test_for_zero(t4ddg_1(2,1,2,0) - (t4ddg_2(2,1,2,0)/t2s_2(2,0))
		,"T2s(k,l)%T4ddg(i,j,k,l)(2,1,2,0)");
  test_for_zero(t4ddg_1(2,1,2,1) - (t4ddg_2(2,1,2,1)/t2s_2(2,1))
		,"T2s(k,l)%T4ddg(i,j,k,l)(2,1,2,1)");
  test_for_zero(t4ddg_1(2,1,2,2) - (t4ddg_2(2,1,2,2)/t2s_2(2,2))
		,"T2s(k,l)%T4ddg(i,j,k,l)(2,1,2,2)");
  test_for_zero(t4ddg_1(2,2,0,0) - (t4ddg_2(2,2,0,0)/t2s_2(0,0))
		,"T2s(k,l)%T4ddg(i,j,k,l)(2,2,0,0)");
  test_for_zero(t4ddg_1(2,2,0,1) - (t4ddg_2(2,2,0,1)/t2s_2(0,1))
		,"T2s(k,l)%T4ddg(i,j,k,l)(2,2,0,1)");
  test_for_zero(t4ddg_1(2,2,0,2) - (t4ddg_2(2,2,0,2)/t2s_2(0,2))
		,"T2s(k,l)%T4ddg(i,j,k,l)(2,2,0,2)");
  test_for_zero(t4ddg_1(2,2,1,0) - (t4ddg_2(2,2,1,0)/t2s_2(1,0))
		,"T2s(k,l)%T4ddg(i,j,k,l)(2,2,1,0)");
  test_for_zero(t4ddg_1(2,2,1,1) - (t4ddg_2(2,2,1,1)/t2s_2(1,1))
		,"T2s(k,l)%T4ddg(i,j,k,l)(2,2,1,1)");
  test_for_zero(t4ddg_1(2,2,1,2) - (t4ddg_2(2,2,1,2)/t2s_2(1,2))
		,"T2s(k,l)%T4ddg(i,j,k,l)(2,2,1,2)");
  test_for_zero(t4ddg_1(2,2,2,0) - (t4ddg_2(2,2,2,0)/t2s_2(2,0))
		,"T2s(k,l)%T4ddg(i,j,k,l)(2,2,2,0)");
  test_for_zero(t4ddg_1(2,2,2,1) - (t4ddg_2(2,2,2,1)/t2s_2(2,1))
		,"T2s(k,l)%T4ddg(i,j,k,l)(2,2,2,1)");
  test_for_zero(t4ddg_1(2,2,2,2) - (t4ddg_2(2,2,2,2)/t2s_2(2,2))
		,"T2s(k,l)%T4ddg(i,j,k,l)(2,2,2,2)");

  /* I originally put these declarations for unknown reasons, but they
     won't work because the result is not a Tensor4_ddg.  The
     multiplication messes up the symmetries. */


//    t4ddg_1(i,j,k,l)=(t4ddg_2(i,j,k,l)&t2s_3(j,l));
//    test_for_zero(t4ddg_1(0,0,0,0) - (t4ddg_2(0,0,0,0)*t2s_3(0,0))
//        test_for_zero(t4ddg_1(0,0,0,1) - (t4ddg_2(0,0,0,1)*t2s_3(0,1))
//        test_for_zero(t4ddg_1(0,0,0,2) - (t4ddg_2(0,0,0,2)*t2s_3(0,2))
//        test_for_zero(t4ddg_1(0,0,1,0) - (t4ddg_2(0,0,1,0)*t2s_3(0,0))
//        test_for_zero(t4ddg_1(0,0,1,1) - (t4ddg_2(0,0,1,1)*t2s_3(0,1))
//        test_for_zero(t4ddg_1(0,0,1,2) - (t4ddg_2(0,0,1,2)*t2s_3(0,2))
//        test_for_zero(t4ddg_1(0,0,2,0) - (t4ddg_2(0,0,2,0)*t2s_3(0,0))
//        test_for_zero(t4ddg_1(0,0,2,1) - (t4ddg_2(0,0,2,1)*t2s_3(0,1))
//        test_for_zero(t4ddg_1(0,0,2,2) - (t4ddg_2(0,0,2,2)*t2s_3(0,2))
//        test_for_zero(t4ddg_1(0,1,0,0) - (t4ddg_2(0,1,0,0)*t2s_3(1,0))
//        test_for_zero(t4ddg_1(0,1,0,1) - (t4ddg_2(0,1,0,1)*t2s_3(1,1))
//        test_for_zero(t4ddg_1(0,1,0,2) - (t4ddg_2(0,1,0,2)*t2s_3(1,2))
//        test_for_zero(t4ddg_1(0,1,1,0) - (t4ddg_2(0,1,1,0)*t2s_3(1,0))
//        test_for_zero(t4ddg_1(0,1,1,1) - (t4ddg_2(0,1,1,1)*t2s_3(1,1))
//        test_for_zero(t4ddg_1(0,1,1,2) - (t4ddg_2(0,1,1,2)*t2s_3(1,2))
//        test_for_zero(t4ddg_1(0,1,2,0) - (t4ddg_2(0,1,2,0)*t2s_3(1,0))
//        test_for_zero(t4ddg_1(0,1,2,1) - (t4ddg_2(0,1,2,1)*t2s_3(1,1))
//        test_for_zero(t4ddg_1(0,1,2,2) - (t4ddg_2(0,1,2,2)*t2s_3(1,2))
//        test_for_zero(t4ddg_1(0,2,0,0) - (t4ddg_2(0,2,0,0)*t2s_3(2,0))
//        test_for_zero(t4ddg_1(0,2,0,1) - (t4ddg_2(0,2,0,1)*t2s_3(2,1))
//        test_for_zero(t4ddg_1(0,2,0,2) - (t4ddg_2(0,2,0,2)*t2s_3(2,2))
//        test_for_zero(t4ddg_1(0,2,1,0) - (t4ddg_2(0,2,1,0)*t2s_3(2,0))
//        test_for_zero(t4ddg_1(0,2,1,1) - (t4ddg_2(0,2,1,1)*t2s_3(2,1))
//        test_for_zero(t4ddg_1(0,2,1,2) - (t4ddg_2(0,2,1,2)*t2s_3(2,2))
//        test_for_zero(t4ddg_1(0,2,2,0) - (t4ddg_2(0,2,2,0)*t2s_3(2,0))
//        test_for_zero(t4ddg_1(0,2,2,1) - (t4ddg_2(0,2,2,1)*t2s_3(2,1))
//        test_for_zero(t4ddg_1(0,2,2,2) - (t4ddg_2(0,2,2,2)*t2s_3(2,2))
//        test_for_zero(t4ddg_1(1,0,0,0) - (t4ddg_2(1,0,0,0)*t2s_3(0,0))
//        test_for_zero(t4ddg_1(1,0,0,1) - (t4ddg_2(1,0,0,1)*t2s_3(0,1))
//        test_for_zero(t4ddg_1(1,0,0,2) - (t4ddg_2(1,0,0,2)*t2s_3(0,2))
//        test_for_zero(t4ddg_1(1,0,1,0) - (t4ddg_2(1,0,1,0)*t2s_3(0,0))
//        test_for_zero(t4ddg_1(1,0,1,1) - (t4ddg_2(1,0,1,1)*t2s_3(0,1))
//        test_for_zero(t4ddg_1(1,0,1,2) - (t4ddg_2(1,0,1,2)*t2s_3(0,2))
//        test_for_zero(t4ddg_1(1,0,2,0) - (t4ddg_2(1,0,2,0)*t2s_3(0,0))
//        test_for_zero(t4ddg_1(1,0,2,1) - (t4ddg_2(1,0,2,1)*t2s_3(0,1))
//        test_for_zero(t4ddg_1(1,0,2,2) - (t4ddg_2(1,0,2,2)*t2s_3(0,2))
//        test_for_zero(t4ddg_1(1,1,0,0) - (t4ddg_2(1,1,0,0)*t2s_3(1,0))
//        test_for_zero(t4ddg_1(1,1,0,1) - (t4ddg_2(1,1,0,1)*t2s_3(1,1))
//        test_for_zero(t4ddg_1(1,1,0,2) - (t4ddg_2(1,1,0,2)*t2s_3(1,2))
//        test_for_zero(t4ddg_1(1,1,1,0) - (t4ddg_2(1,1,1,0)*t2s_3(1,0))
//        test_for_zero(t4ddg_1(1,1,1,1) - (t4ddg_2(1,1,1,1)*t2s_3(1,1))
//        test_for_zero(t4ddg_1(1,1,1,2) - (t4ddg_2(1,1,1,2)*t2s_3(1,2))
//        test_for_zero(t4ddg_1(1,1,2,0) - (t4ddg_2(1,1,2,0)*t2s_3(1,0))
//        test_for_zero(t4ddg_1(1,1,2,1) - (t4ddg_2(1,1,2,1)*t2s_3(1,1))
//        test_for_zero(t4ddg_1(1,1,2,2) - (t4ddg_2(1,1,2,2)*t2s_3(1,2))
//        test_for_zero(t4ddg_1(1,2,0,0) - (t4ddg_2(1,2,0,0)*t2s_3(2,0))
//        test_for_zero(t4ddg_1(1,2,0,1) - (t4ddg_2(1,2,0,1)*t2s_3(2,1))
//        test_for_zero(t4ddg_1(1,2,0,2) - (t4ddg_2(1,2,0,2)*t2s_3(2,2))
//        test_for_zero(t4ddg_1(1,2,1,0) - (t4ddg_2(1,2,1,0)*t2s_3(2,0))
//        test_for_zero(t4ddg_1(1,2,1,1) - (t4ddg_2(1,2,1,1)*t2s_3(2,1))
//        test_for_zero(t4ddg_1(1,2,1,2) - (t4ddg_2(1,2,1,2)*t2s_3(2,2))
//        test_for_zero(t4ddg_1(1,2,2,0) - (t4ddg_2(1,2,2,0)*t2s_3(2,0))
//        test_for_zero(t4ddg_1(1,2,2,1) - (t4ddg_2(1,2,2,1)*t2s_3(2,1))
//        test_for_zero(t4ddg_1(1,2,2,2) - (t4ddg_2(1,2,2,2)*t2s_3(2,2))
//        test_for_zero(t4ddg_1(2,0,0,0) - (t4ddg_2(2,0,0,0)*t2s_3(0,0))
//        test_for_zero(t4ddg_1(2,0,0,1) - (t4ddg_2(2,0,0,1)*t2s_3(0,1))
//        test_for_zero(t4ddg_1(2,0,0,2) - (t4ddg_2(2,0,0,2)*t2s_3(0,2))
//        test_for_zero(t4ddg_1(2,0,1,0) - (t4ddg_2(2,0,1,0)*t2s_3(0,0))
//        test_for_zero(t4ddg_1(2,0,1,1) - (t4ddg_2(2,0,1,1)*t2s_3(0,1))
//        test_for_zero(t4ddg_1(2,0,1,2) - (t4ddg_2(2,0,1,2)*t2s_3(0,2))
//        test_for_zero(t4ddg_1(2,0,2,0) - (t4ddg_2(2,0,2,0)*t2s_3(0,0))
//        test_for_zero(t4ddg_1(2,0,2,1) - (t4ddg_2(2,0,2,1)*t2s_3(0,1))
//        test_for_zero(t4ddg_1(2,0,2,2) - (t4ddg_2(2,0,2,2)*t2s_3(0,2))
//        test_for_zero(t4ddg_1(2,1,0,0) - (t4ddg_2(2,1,0,0)*t2s_3(1,0))
//        test_for_zero(t4ddg_1(2,1,0,1) - (t4ddg_2(2,1,0,1)*t2s_3(1,1))
//        test_for_zero(t4ddg_1(2,1,0,2) - (t4ddg_2(2,1,0,2)*t2s_3(1,2))
//        test_for_zero(t4ddg_1(2,1,1,0) - (t4ddg_2(2,1,1,0)*t2s_3(1,0))
//        test_for_zero(t4ddg_1(2,1,1,1) - (t4ddg_2(2,1,1,1)*t2s_3(1,1))
//        test_for_zero(t4ddg_1(2,1,1,2) - (t4ddg_2(2,1,1,2)*t2s_3(1,2))
//        test_for_zero(t4ddg_1(2,1,2,0) - (t4ddg_2(2,1,2,0)*t2s_3(1,0))
//        test_for_zero(t4ddg_1(2,1,2,1) - (t4ddg_2(2,1,2,1)*t2s_3(1,1))
//        test_for_zero(t4ddg_1(2,1,2,2) - (t4ddg_2(2,1,2,2)*t2s_3(1,2))
//        test_for_zero(t4ddg_1(2,2,0,0) - (t4ddg_2(2,2,0,0)*t2s_3(2,0))
//        test_for_zero(t4ddg_1(2,2,0,1) - (t4ddg_2(2,2,0,1)*t2s_3(2,1))
//        test_for_zero(t4ddg_1(2,2,0,2) - (t4ddg_2(2,2,0,2)*t2s_3(2,2))
//        test_for_zero(t4ddg_1(2,2,1,0) - (t4ddg_2(2,2,1,0)*t2s_3(2,0))
//        test_for_zero(t4ddg_1(2,2,1,1) - (t4ddg_2(2,2,1,1)*t2s_3(2,1))
//        test_for_zero(t4ddg_1(2,2,1,2) - (t4ddg_2(2,2,1,2)*t2s_3(2,2))
//        test_for_zero(t4ddg_1(2,2,2,0) - (t4ddg_2(2,2,2,0)*t2s_3(2,0))
//        test_for_zero(t4ddg_1(2,2,2,1) - (t4ddg_2(2,2,2,1)*t2s_3(2,1))
//        test_for_zero(t4ddg_1(2,2,2,2) - (t4ddg_2(2,2,2,2)*t2s_3(2,2)));


//    t4ddg_1(i,j,k,l)=(t2s_2(j,l)&t4ddg_2(i,j,k,l));
//    test_for_zero(t4ddg_1(0,0,0,0) - (t4ddg_2(0,0,0,0)*t2s_2(0,0))
//        test_for_zero(t4ddg_1(0,0,0,1) - (t4ddg_2(0,0,0,1)*t2s_2(0,1))
//        test_for_zero(t4ddg_1(0,0,0,2) - (t4ddg_2(0,0,0,2)*t2s_2(0,2))
//        test_for_zero(t4ddg_1(0,0,1,0) - (t4ddg_2(0,0,1,0)*t2s_2(0,0))
//        test_for_zero(t4ddg_1(0,0,1,1) - (t4ddg_2(0,0,1,1)*t2s_2(0,1))
//        test_for_zero(t4ddg_1(0,0,1,2) - (t4ddg_2(0,0,1,2)*t2s_2(0,2))
//        test_for_zero(t4ddg_1(0,0,2,0) - (t4ddg_2(0,0,2,0)*t2s_2(0,0))
//        test_for_zero(t4ddg_1(0,0,2,1) - (t4ddg_2(0,0,2,1)*t2s_2(0,1))
//        test_for_zero(t4ddg_1(0,0,2,2) - (t4ddg_2(0,0,2,2)*t2s_2(0,2))
//        test_for_zero(t4ddg_1(0,1,0,0) - (t4ddg_2(0,1,0,0)*t2s_2(1,0))
//        test_for_zero(t4ddg_1(0,1,0,1) - (t4ddg_2(0,1,0,1)*t2s_2(1,1))
//        test_for_zero(t4ddg_1(0,1,0,2) - (t4ddg_2(0,1,0,2)*t2s_2(1,2))
//        test_for_zero(t4ddg_1(0,1,1,0) - (t4ddg_2(0,1,1,0)*t2s_2(1,0))
//        test_for_zero(t4ddg_1(0,1,1,1) - (t4ddg_2(0,1,1,1)*t2s_2(1,1))
//        test_for_zero(t4ddg_1(0,1,1,2) - (t4ddg_2(0,1,1,2)*t2s_2(1,2))
//        test_for_zero(t4ddg_1(0,1,2,0) - (t4ddg_2(0,1,2,0)*t2s_2(1,0))
//        test_for_zero(t4ddg_1(0,1,2,1) - (t4ddg_2(0,1,2,1)*t2s_2(1,1))
//        test_for_zero(t4ddg_1(0,1,2,2) - (t4ddg_2(0,1,2,2)*t2s_2(1,2))
//        test_for_zero(t4ddg_1(0,2,0,0) - (t4ddg_2(0,2,0,0)*t2s_2(2,0))
//        test_for_zero(t4ddg_1(0,2,0,1) - (t4ddg_2(0,2,0,1)*t2s_2(2,1))
//        test_for_zero(t4ddg_1(0,2,0,2) - (t4ddg_2(0,2,0,2)*t2s_2(2,2))
//        test_for_zero(t4ddg_1(0,2,1,0) - (t4ddg_2(0,2,1,0)*t2s_2(2,0))
//        test_for_zero(t4ddg_1(0,2,1,1) - (t4ddg_2(0,2,1,1)*t2s_2(2,1))
//        test_for_zero(t4ddg_1(0,2,1,2) - (t4ddg_2(0,2,1,2)*t2s_2(2,2))
//        test_for_zero(t4ddg_1(0,2,2,0) - (t4ddg_2(0,2,2,0)*t2s_2(2,0))
//        test_for_zero(t4ddg_1(0,2,2,1) - (t4ddg_2(0,2,2,1)*t2s_2(2,1))
//        test_for_zero(t4ddg_1(0,2,2,2) - (t4ddg_2(0,2,2,2)*t2s_2(2,2))
//        test_for_zero(t4ddg_1(1,0,0,0) - (t4ddg_2(1,0,0,0)*t2s_2(0,0))
//        test_for_zero(t4ddg_1(1,0,0,1) - (t4ddg_2(1,0,0,1)*t2s_2(0,1))
//        test_for_zero(t4ddg_1(1,0,0,2) - (t4ddg_2(1,0,0,2)*t2s_2(0,2))
//        test_for_zero(t4ddg_1(1,0,1,0) - (t4ddg_2(1,0,1,0)*t2s_2(0,0))
//        test_for_zero(t4ddg_1(1,0,1,1) - (t4ddg_2(1,0,1,1)*t2s_2(0,1))
//        test_for_zero(t4ddg_1(1,0,1,2) - (t4ddg_2(1,0,1,2)*t2s_2(0,2))
//        test_for_zero(t4ddg_1(1,0,2,0) - (t4ddg_2(1,0,2,0)*t2s_2(0,0))
//        test_for_zero(t4ddg_1(1,0,2,1) - (t4ddg_2(1,0,2,1)*t2s_2(0,1))
//        test_for_zero(t4ddg_1(1,0,2,2) - (t4ddg_2(1,0,2,2)*t2s_2(0,2))
//        test_for_zero(t4ddg_1(1,1,0,0) - (t4ddg_2(1,1,0,0)*t2s_2(1,0))
//        test_for_zero(t4ddg_1(1,1,0,1) - (t4ddg_2(1,1,0,1)*t2s_2(1,1))
//        test_for_zero(t4ddg_1(1,1,0,2) - (t4ddg_2(1,1,0,2)*t2s_2(1,2))
//        test_for_zero(t4ddg_1(1,1,1,0) - (t4ddg_2(1,1,1,0)*t2s_2(1,0))
//        test_for_zero(t4ddg_1(1,1,1,1) - (t4ddg_2(1,1,1,1)*t2s_2(1,1))
//        test_for_zero(t4ddg_1(1,1,1,2) - (t4ddg_2(1,1,1,2)*t2s_2(1,2))
//        test_for_zero(t4ddg_1(1,1,2,0) - (t4ddg_2(1,1,2,0)*t2s_2(1,0))
//        test_for_zero(t4ddg_1(1,1,2,1) - (t4ddg_2(1,1,2,1)*t2s_2(1,1))
//        test_for_zero(t4ddg_1(1,1,2,2) - (t4ddg_2(1,1,2,2)*t2s_2(1,2))
//        test_for_zero(t4ddg_1(1,2,0,0) - (t4ddg_2(1,2,0,0)*t2s_2(2,0))
//        test_for_zero(t4ddg_1(1,2,0,1) - (t4ddg_2(1,2,0,1)*t2s_2(2,1))
//        test_for_zero(t4ddg_1(1,2,0,2) - (t4ddg_2(1,2,0,2)*t2s_2(2,2))
//        test_for_zero(t4ddg_1(1,2,1,0) - (t4ddg_2(1,2,1,0)*t2s_2(2,0))
//        test_for_zero(t4ddg_1(1,2,1,1) - (t4ddg_2(1,2,1,1)*t2s_2(2,1))
//        test_for_zero(t4ddg_1(1,2,1,2) - (t4ddg_2(1,2,1,2)*t2s_2(2,2))
//        test_for_zero(t4ddg_1(1,2,2,0) - (t4ddg_2(1,2,2,0)*t2s_2(2,0))
//        test_for_zero(t4ddg_1(1,2,2,1) - (t4ddg_2(1,2,2,1)*t2s_2(2,1))
//        test_for_zero(t4ddg_1(1,2,2,2) - (t4ddg_2(1,2,2,2)*t2s_2(2,2))
//        test_for_zero(t4ddg_1(2,0,0,0) - (t4ddg_2(2,0,0,0)*t2s_2(0,0))
//        test_for_zero(t4ddg_1(2,0,0,1) - (t4ddg_2(2,0,0,1)*t2s_2(0,1))
//        test_for_zero(t4ddg_1(2,0,0,2) - (t4ddg_2(2,0,0,2)*t2s_2(0,2))
//        test_for_zero(t4ddg_1(2,0,1,0) - (t4ddg_2(2,0,1,0)*t2s_2(0,0))
//        test_for_zero(t4ddg_1(2,0,1,1) - (t4ddg_2(2,0,1,1)*t2s_2(0,1))
//        test_for_zero(t4ddg_1(2,0,1,2) - (t4ddg_2(2,0,1,2)*t2s_2(0,2))
//        test_for_zero(t4ddg_1(2,0,2,0) - (t4ddg_2(2,0,2,0)*t2s_2(0,0))
//        test_for_zero(t4ddg_1(2,0,2,1) - (t4ddg_2(2,0,2,1)*t2s_2(0,1))
//        test_for_zero(t4ddg_1(2,0,2,2) - (t4ddg_2(2,0,2,2)*t2s_2(0,2))
//        test_for_zero(t4ddg_1(2,1,0,0) - (t4ddg_2(2,1,0,0)*t2s_2(1,0))
//        test_for_zero(t4ddg_1(2,1,0,1) - (t4ddg_2(2,1,0,1)*t2s_2(1,1))
//        test_for_zero(t4ddg_1(2,1,0,2) - (t4ddg_2(2,1,0,2)*t2s_2(1,2))
//        test_for_zero(t4ddg_1(2,1,1,0) - (t4ddg_2(2,1,1,0)*t2s_2(1,0))
//        test_for_zero(t4ddg_1(2,1,1,1) - (t4ddg_2(2,1,1,1)*t2s_2(1,1))
//        test_for_zero(t4ddg_1(2,1,1,2) - (t4ddg_2(2,1,1,2)*t2s_2(1,2))
//        test_for_zero(t4ddg_1(2,1,2,0) - (t4ddg_2(2,1,2,0)*t2s_2(1,0))
//        test_for_zero(t4ddg_1(2,1,2,1) - (t4ddg_2(2,1,2,1)*t2s_2(1,1))
//        test_for_zero(t4ddg_1(2,1,2,2) - (t4ddg_2(2,1,2,2)*t2s_2(1,2))
//        test_for_zero(t4ddg_1(2,2,0,0) - (t4ddg_2(2,2,0,0)*t2s_2(2,0))
//        test_for_zero(t4ddg_1(2,2,0,1) - (t4ddg_2(2,2,0,1)*t2s_2(2,1))
//        test_for_zero(t4ddg_1(2,2,0,2) - (t4ddg_2(2,2,0,2)*t2s_2(2,2))
//        test_for_zero(t4ddg_1(2,2,1,0) - (t4ddg_2(2,2,1,0)*t2s_2(2,0))
//        test_for_zero(t4ddg_1(2,2,1,1) - (t4ddg_2(2,2,1,1)*t2s_2(2,1))
//        test_for_zero(t4ddg_1(2,2,1,2) - (t4ddg_2(2,2,1,2)*t2s_2(2,2))
//        test_for_zero(t4ddg_1(2,2,2,0) - (t4ddg_2(2,2,2,0)*t2s_2(2,0))
//        test_for_zero(t4ddg_1(2,2,2,1) - (t4ddg_2(2,2,2,1)*t2s_2(2,1))
//        test_for_zero(t4ddg_1(2,2,2,2) - (t4ddg_2(2,2,2,2)*t2s_2(2,2)));
//    cout << endl;


//    t4ddg_1(i,j,k,l)=(t4ddg_2(i,j,k,l)&t2s_3(l,j));
//    test_for_zero(t4ddg_1(0,0,0,0) - (t4ddg_2(0,0,0,0)*t2s_3(0,0))
//        test_for_zero(t4ddg_1(0,0,0,1) - (t4ddg_2(0,0,0,1)*t2s_3(0,1))
//        test_for_zero(t4ddg_1(0,0,0,2) - (t4ddg_2(0,0,0,2)*t2s_3(0,2))
//        test_for_zero(t4ddg_1(0,0,1,0) - (t4ddg_2(0,0,1,0)*t2s_3(0,0))
//        test_for_zero(t4ddg_1(0,0,1,1) - (t4ddg_2(0,0,1,1)*t2s_3(0,1))
//        test_for_zero(t4ddg_1(0,0,1,2) - (t4ddg_2(0,0,1,2)*t2s_3(0,2))
//        test_for_zero(t4ddg_1(0,0,2,0) - (t4ddg_2(0,0,2,0)*t2s_3(0,0))
//        test_for_zero(t4ddg_1(0,0,2,1) - (t4ddg_2(0,0,2,1)*t2s_3(0,1))
//        test_for_zero(t4ddg_1(0,0,2,2) - (t4ddg_2(0,0,2,2)*t2s_3(0,2))
//        test_for_zero(t4ddg_1(0,1,0,0) - (t4ddg_2(0,1,0,0)*t2s_3(1,0))
//        test_for_zero(t4ddg_1(0,1,0,1) - (t4ddg_2(0,1,0,1)*t2s_3(1,1))
//        test_for_zero(t4ddg_1(0,1,0,2) - (t4ddg_2(0,1,0,2)*t2s_3(1,2))
//        test_for_zero(t4ddg_1(0,1,1,0) - (t4ddg_2(0,1,1,0)*t2s_3(1,0))
//        test_for_zero(t4ddg_1(0,1,1,1) - (t4ddg_2(0,1,1,1)*t2s_3(1,1))
//        test_for_zero(t4ddg_1(0,1,1,2) - (t4ddg_2(0,1,1,2)*t2s_3(1,2))
//        test_for_zero(t4ddg_1(0,1,2,0) - (t4ddg_2(0,1,2,0)*t2s_3(1,0))
//        test_for_zero(t4ddg_1(0,1,2,1) - (t4ddg_2(0,1,2,1)*t2s_3(1,1))
//        test_for_zero(t4ddg_1(0,1,2,2) - (t4ddg_2(0,1,2,2)*t2s_3(1,2))
//        test_for_zero(t4ddg_1(0,2,0,0) - (t4ddg_2(0,2,0,0)*t2s_3(2,0))
//        test_for_zero(t4ddg_1(0,2,0,1) - (t4ddg_2(0,2,0,1)*t2s_3(2,1))
//        test_for_zero(t4ddg_1(0,2,0,2) - (t4ddg_2(0,2,0,2)*t2s_3(2,2))
//        test_for_zero(t4ddg_1(0,2,1,0) - (t4ddg_2(0,2,1,0)*t2s_3(2,0))
//        test_for_zero(t4ddg_1(0,2,1,1) - (t4ddg_2(0,2,1,1)*t2s_3(2,1))
//        test_for_zero(t4ddg_1(0,2,1,2) - (t4ddg_2(0,2,1,2)*t2s_3(2,2))
//        test_for_zero(t4ddg_1(0,2,2,0) - (t4ddg_2(0,2,2,0)*t2s_3(2,0))
//        test_for_zero(t4ddg_1(0,2,2,1) - (t4ddg_2(0,2,2,1)*t2s_3(2,1))
//        test_for_zero(t4ddg_1(0,2,2,2) - (t4ddg_2(0,2,2,2)*t2s_3(2,2))
//        test_for_zero(t4ddg_1(1,0,0,0) - (t4ddg_2(1,0,0,0)*t2s_3(0,0))
//        test_for_zero(t4ddg_1(1,0,0,1) - (t4ddg_2(1,0,0,1)*t2s_3(0,1))
//        test_for_zero(t4ddg_1(1,0,0,2) - (t4ddg_2(1,0,0,2)*t2s_3(0,2))
//        test_for_zero(t4ddg_1(1,0,1,0) - (t4ddg_2(1,0,1,0)*t2s_3(0,0))
//        test_for_zero(t4ddg_1(1,0,1,1) - (t4ddg_2(1,0,1,1)*t2s_3(0,1))
//        test_for_zero(t4ddg_1(1,0,1,2) - (t4ddg_2(1,0,1,2)*t2s_3(0,2))
//        test_for_zero(t4ddg_1(1,0,2,0) - (t4ddg_2(1,0,2,0)*t2s_3(0,0))
//        test_for_zero(t4ddg_1(1,0,2,1) - (t4ddg_2(1,0,2,1)*t2s_3(0,1))
//        test_for_zero(t4ddg_1(1,0,2,2) - (t4ddg_2(1,0,2,2)*t2s_3(0,2))
//        test_for_zero(t4ddg_1(1,1,0,0) - (t4ddg_2(1,1,0,0)*t2s_3(1,0))
//        test_for_zero(t4ddg_1(1,1,0,1) - (t4ddg_2(1,1,0,1)*t2s_3(1,1))
//        test_for_zero(t4ddg_1(1,1,0,2) - (t4ddg_2(1,1,0,2)*t2s_3(1,2))
//        test_for_zero(t4ddg_1(1,1,1,0) - (t4ddg_2(1,1,1,0)*t2s_3(1,0))
//        test_for_zero(t4ddg_1(1,1,1,1) - (t4ddg_2(1,1,1,1)*t2s_3(1,1))
//        test_for_zero(t4ddg_1(1,1,1,2) - (t4ddg_2(1,1,1,2)*t2s_3(1,2))
//        test_for_zero(t4ddg_1(1,1,2,0) - (t4ddg_2(1,1,2,0)*t2s_3(1,0))
//        test_for_zero(t4ddg_1(1,1,2,1) - (t4ddg_2(1,1,2,1)*t2s_3(1,1))
//        test_for_zero(t4ddg_1(1,1,2,2) - (t4ddg_2(1,1,2,2)*t2s_3(1,2))
//        test_for_zero(t4ddg_1(1,2,0,0) - (t4ddg_2(1,2,0,0)*t2s_3(2,0))
//        test_for_zero(t4ddg_1(1,2,0,1) - (t4ddg_2(1,2,0,1)*t2s_3(2,1))
//        test_for_zero(t4ddg_1(1,2,0,2) - (t4ddg_2(1,2,0,2)*t2s_3(2,2))
//        test_for_zero(t4ddg_1(1,2,1,0) - (t4ddg_2(1,2,1,0)*t2s_3(2,0))
//        test_for_zero(t4ddg_1(1,2,1,1) - (t4ddg_2(1,2,1,1)*t2s_3(2,1))
//        test_for_zero(t4ddg_1(1,2,1,2) - (t4ddg_2(1,2,1,2)*t2s_3(2,2))
//        test_for_zero(t4ddg_1(1,2,2,0) - (t4ddg_2(1,2,2,0)*t2s_3(2,0))
//        test_for_zero(t4ddg_1(1,2,2,1) - (t4ddg_2(1,2,2,1)*t2s_3(2,1))
//        test_for_zero(t4ddg_1(1,2,2,2) - (t4ddg_2(1,2,2,2)*t2s_3(2,2))
//        test_for_zero(t4ddg_1(2,0,0,0) - (t4ddg_2(2,0,0,0)*t2s_3(0,0))
//        test_for_zero(t4ddg_1(2,0,0,1) - (t4ddg_2(2,0,0,1)*t2s_3(0,1))
//        test_for_zero(t4ddg_1(2,0,0,2) - (t4ddg_2(2,0,0,2)*t2s_3(0,2))
//        test_for_zero(t4ddg_1(2,0,1,0) - (t4ddg_2(2,0,1,0)*t2s_3(0,0))
//        test_for_zero(t4ddg_1(2,0,1,1) - (t4ddg_2(2,0,1,1)*t2s_3(0,1))
//        test_for_zero(t4ddg_1(2,0,1,2) - (t4ddg_2(2,0,1,2)*t2s_3(0,2))
//        test_for_zero(t4ddg_1(2,0,2,0) - (t4ddg_2(2,0,2,0)*t2s_3(0,0))
//        test_for_zero(t4ddg_1(2,0,2,1) - (t4ddg_2(2,0,2,1)*t2s_3(0,1))
//        test_for_zero(t4ddg_1(2,0,2,2) - (t4ddg_2(2,0,2,2)*t2s_3(0,2))
//        test_for_zero(t4ddg_1(2,1,0,0) - (t4ddg_2(2,1,0,0)*t2s_3(1,0))
//        test_for_zero(t4ddg_1(2,1,0,1) - (t4ddg_2(2,1,0,1)*t2s_3(1,1))
//        test_for_zero(t4ddg_1(2,1,0,2) - (t4ddg_2(2,1,0,2)*t2s_3(1,2))
//        test_for_zero(t4ddg_1(2,1,1,0) - (t4ddg_2(2,1,1,0)*t2s_3(1,0))
//        test_for_zero(t4ddg_1(2,1,1,1) - (t4ddg_2(2,1,1,1)*t2s_3(1,1))
//        test_for_zero(t4ddg_1(2,1,1,2) - (t4ddg_2(2,1,1,2)*t2s_3(1,2))
//        test_for_zero(t4ddg_1(2,1,2,0) - (t4ddg_2(2,1,2,0)*t2s_3(1,0))
//        test_for_zero(t4ddg_1(2,1,2,1) - (t4ddg_2(2,1,2,1)*t2s_3(1,1))
//        test_for_zero(t4ddg_1(2,1,2,2) - (t4ddg_2(2,1,2,2)*t2s_3(1,2))
//        test_for_zero(t4ddg_1(2,2,0,0) - (t4ddg_2(2,2,0,0)*t2s_3(2,0))
//        test_for_zero(t4ddg_1(2,2,0,1) - (t4ddg_2(2,2,0,1)*t2s_3(2,1))
//        test_for_zero(t4ddg_1(2,2,0,2) - (t4ddg_2(2,2,0,2)*t2s_3(2,2))
//        test_for_zero(t4ddg_1(2,2,1,0) - (t4ddg_2(2,2,1,0)*t2s_3(2,0))
//        test_for_zero(t4ddg_1(2,2,1,1) - (t4ddg_2(2,2,1,1)*t2s_3(2,1))
//        test_for_zero(t4ddg_1(2,2,1,2) - (t4ddg_2(2,2,1,2)*t2s_3(2,2))
//        test_for_zero(t4ddg_1(2,2,2,0) - (t4ddg_2(2,2,2,0)*t2s_3(2,0))
//        test_for_zero(t4ddg_1(2,2,2,1) - (t4ddg_2(2,2,2,1)*t2s_3(2,1))
//        test_for_zero(t4ddg_1(2,2,2,2) - (t4ddg_2(2,2,2,2)*t2s_3(2,2)));
//    cout << endl;


//    t4ddg_1(i,j,k,l)=(t2s_2(l,j)&t4ddg_2(i,j,k,l));
//    test_for_zero(t4ddg_1(0,0,0,0) - (t4ddg_2(0,0,0,0)*t2s_2(0,0))
//        test_for_zero(t4ddg_1(0,0,0,1) - (t4ddg_2(0,0,0,1)*t2s_2(0,1))
//        test_for_zero(t4ddg_1(0,0,0,2) - (t4ddg_2(0,0,0,2)*t2s_2(0,2))
//        test_for_zero(t4ddg_1(0,0,1,0) - (t4ddg_2(0,0,1,0)*t2s_2(0,0))
//        test_for_zero(t4ddg_1(0,0,1,1) - (t4ddg_2(0,0,1,1)*t2s_2(0,1))
//        test_for_zero(t4ddg_1(0,0,1,2) - (t4ddg_2(0,0,1,2)*t2s_2(0,2))
//        test_for_zero(t4ddg_1(0,0,2,0) - (t4ddg_2(0,0,2,0)*t2s_2(0,0))
//        test_for_zero(t4ddg_1(0,0,2,1) - (t4ddg_2(0,0,2,1)*t2s_2(0,1))
//        test_for_zero(t4ddg_1(0,0,2,2) - (t4ddg_2(0,0,2,2)*t2s_2(0,2))
//        test_for_zero(t4ddg_1(0,1,0,0) - (t4ddg_2(0,1,0,0)*t2s_2(1,0))
//        test_for_zero(t4ddg_1(0,1,0,1) - (t4ddg_2(0,1,0,1)*t2s_2(1,1))
//        test_for_zero(t4ddg_1(0,1,0,2) - (t4ddg_2(0,1,0,2)*t2s_2(1,2))
//        test_for_zero(t4ddg_1(0,1,1,0) - (t4ddg_2(0,1,1,0)*t2s_2(1,0))
//        test_for_zero(t4ddg_1(0,1,1,1) - (t4ddg_2(0,1,1,1)*t2s_2(1,1))
//        test_for_zero(t4ddg_1(0,1,1,2) - (t4ddg_2(0,1,1,2)*t2s_2(1,2))
//        test_for_zero(t4ddg_1(0,1,2,0) - (t4ddg_2(0,1,2,0)*t2s_2(1,0))
//        test_for_zero(t4ddg_1(0,1,2,1) - (t4ddg_2(0,1,2,1)*t2s_2(1,1))
//        test_for_zero(t4ddg_1(0,1,2,2) - (t4ddg_2(0,1,2,2)*t2s_2(1,2))
//        test_for_zero(t4ddg_1(0,2,0,0) - (t4ddg_2(0,2,0,0)*t2s_2(2,0))
//        test_for_zero(t4ddg_1(0,2,0,1) - (t4ddg_2(0,2,0,1)*t2s_2(2,1))
//        test_for_zero(t4ddg_1(0,2,0,2) - (t4ddg_2(0,2,0,2)*t2s_2(2,2))
//        test_for_zero(t4ddg_1(0,2,1,0) - (t4ddg_2(0,2,1,0)*t2s_2(2,0))
//        test_for_zero(t4ddg_1(0,2,1,1) - (t4ddg_2(0,2,1,1)*t2s_2(2,1))
//        test_for_zero(t4ddg_1(0,2,1,2) - (t4ddg_2(0,2,1,2)*t2s_2(2,2))
//        test_for_zero(t4ddg_1(0,2,2,0) - (t4ddg_2(0,2,2,0)*t2s_2(2,0))
//        test_for_zero(t4ddg_1(0,2,2,1) - (t4ddg_2(0,2,2,1)*t2s_2(2,1))
//        test_for_zero(t4ddg_1(0,2,2,2) - (t4ddg_2(0,2,2,2)*t2s_2(2,2))
//        test_for_zero(t4ddg_1(1,0,0,0) - (t4ddg_2(1,0,0,0)*t2s_2(0,0))
//        test_for_zero(t4ddg_1(1,0,0,1) - (t4ddg_2(1,0,0,1)*t2s_2(0,1))
//        test_for_zero(t4ddg_1(1,0,0,2) - (t4ddg_2(1,0,0,2)*t2s_2(0,2))
//        test_for_zero(t4ddg_1(1,0,1,0) - (t4ddg_2(1,0,1,0)*t2s_2(0,0))
//        test_for_zero(t4ddg_1(1,0,1,1) - (t4ddg_2(1,0,1,1)*t2s_2(0,1))
//        test_for_zero(t4ddg_1(1,0,1,2) - (t4ddg_2(1,0,1,2)*t2s_2(0,2))
//        test_for_zero(t4ddg_1(1,0,2,0) - (t4ddg_2(1,0,2,0)*t2s_2(0,0))
//        test_for_zero(t4ddg_1(1,0,2,1) - (t4ddg_2(1,0,2,1)*t2s_2(0,1))
//        test_for_zero(t4ddg_1(1,0,2,2) - (t4ddg_2(1,0,2,2)*t2s_2(0,2))
//        test_for_zero(t4ddg_1(1,1,0,0) - (t4ddg_2(1,1,0,0)*t2s_2(1,0))
//        test_for_zero(t4ddg_1(1,1,0,1) - (t4ddg_2(1,1,0,1)*t2s_2(1,1))
//        test_for_zero(t4ddg_1(1,1,0,2) - (t4ddg_2(1,1,0,2)*t2s_2(1,2))
//        test_for_zero(t4ddg_1(1,1,1,0) - (t4ddg_2(1,1,1,0)*t2s_2(1,0))
//        test_for_zero(t4ddg_1(1,1,1,1) - (t4ddg_2(1,1,1,1)*t2s_2(1,1))
//        test_for_zero(t4ddg_1(1,1,1,2) - (t4ddg_2(1,1,1,2)*t2s_2(1,2))
//        test_for_zero(t4ddg_1(1,1,2,0) - (t4ddg_2(1,1,2,0)*t2s_2(1,0))
//        test_for_zero(t4ddg_1(1,1,2,1) - (t4ddg_2(1,1,2,1)*t2s_2(1,1))
//        test_for_zero(t4ddg_1(1,1,2,2) - (t4ddg_2(1,1,2,2)*t2s_2(1,2))
//        test_for_zero(t4ddg_1(1,2,0,0) - (t4ddg_2(1,2,0,0)*t2s_2(2,0))
//        test_for_zero(t4ddg_1(1,2,0,1) - (t4ddg_2(1,2,0,1)*t2s_2(2,1))
//        test_for_zero(t4ddg_1(1,2,0,2) - (t4ddg_2(1,2,0,2)*t2s_2(2,2))
//        test_for_zero(t4ddg_1(1,2,1,0) - (t4ddg_2(1,2,1,0)*t2s_2(2,0))
//        test_for_zero(t4ddg_1(1,2,1,1) - (t4ddg_2(1,2,1,1)*t2s_2(2,1))
//        test_for_zero(t4ddg_1(1,2,1,2) - (t4ddg_2(1,2,1,2)*t2s_2(2,2))
//        test_for_zero(t4ddg_1(1,2,2,0) - (t4ddg_2(1,2,2,0)*t2s_2(2,0))
//        test_for_zero(t4ddg_1(1,2,2,1) - (t4ddg_2(1,2,2,1)*t2s_2(2,1))
//        test_for_zero(t4ddg_1(1,2,2,2) - (t4ddg_2(1,2,2,2)*t2s_2(2,2))
//        test_for_zero(t4ddg_1(2,0,0,0) - (t4ddg_2(2,0,0,0)*t2s_2(0,0))
//        test_for_zero(t4ddg_1(2,0,0,1) - (t4ddg_2(2,0,0,1)*t2s_2(0,1))
//        test_for_zero(t4ddg_1(2,0,0,2) - (t4ddg_2(2,0,0,2)*t2s_2(0,2))
//        test_for_zero(t4ddg_1(2,0,1,0) - (t4ddg_2(2,0,1,0)*t2s_2(0,0))
//        test_for_zero(t4ddg_1(2,0,1,1) - (t4ddg_2(2,0,1,1)*t2s_2(0,1))
//        test_for_zero(t4ddg_1(2,0,1,2) - (t4ddg_2(2,0,1,2)*t2s_2(0,2))
//        test_for_zero(t4ddg_1(2,0,2,0) - (t4ddg_2(2,0,2,0)*t2s_2(0,0))
//        test_for_zero(t4ddg_1(2,0,2,1) - (t4ddg_2(2,0,2,1)*t2s_2(0,1))
//        test_for_zero(t4ddg_1(2,0,2,2) - (t4ddg_2(2,0,2,2)*t2s_2(0,2))
//        test_for_zero(t4ddg_1(2,1,0,0) - (t4ddg_2(2,1,0,0)*t2s_2(1,0))
//        test_for_zero(t4ddg_1(2,1,0,1) - (t4ddg_2(2,1,0,1)*t2s_2(1,1))
//        test_for_zero(t4ddg_1(2,1,0,2) - (t4ddg_2(2,1,0,2)*t2s_2(1,2))
//        test_for_zero(t4ddg_1(2,1,1,0) - (t4ddg_2(2,1,1,0)*t2s_2(1,0))
//        test_for_zero(t4ddg_1(2,1,1,1) - (t4ddg_2(2,1,1,1)*t2s_2(1,1))
//        test_for_zero(t4ddg_1(2,1,1,2) - (t4ddg_2(2,1,1,2)*t2s_2(1,2))
//        test_for_zero(t4ddg_1(2,1,2,0) - (t4ddg_2(2,1,2,0)*t2s_2(1,0))
//        test_for_zero(t4ddg_1(2,1,2,1) - (t4ddg_2(2,1,2,1)*t2s_2(1,1))
//        test_for_zero(t4ddg_1(2,1,2,2) - (t4ddg_2(2,1,2,2)*t2s_2(1,2))
//        test_for_zero(t4ddg_1(2,2,0,0) - (t4ddg_2(2,2,0,0)*t2s_2(2,0))
//        test_for_zero(t4ddg_1(2,2,0,1) - (t4ddg_2(2,2,0,1)*t2s_2(2,1))
//        test_for_zero(t4ddg_1(2,2,0,2) - (t4ddg_2(2,2,0,2)*t2s_2(2,2))
//        test_for_zero(t4ddg_1(2,2,1,0) - (t4ddg_2(2,2,1,0)*t2s_2(2,0))
//        test_for_zero(t4ddg_1(2,2,1,1) - (t4ddg_2(2,2,1,1)*t2s_2(2,1))
//        test_for_zero(t4ddg_1(2,2,1,2) - (t4ddg_2(2,2,1,2)*t2s_2(2,2))
//        test_for_zero(t4ddg_1(2,2,2,0) - (t4ddg_2(2,2,2,0)*t2s_2(2,0))
//        test_for_zero(t4ddg_1(2,2,2,1) - (t4ddg_2(2,2,2,1)*t2s_2(2,1))
//        test_for_zero(t4ddg_1(2,2,2,2) - (t4ddg_2(2,2,2,2)*t2s_2(2,2)));

  cout << endl;
}
