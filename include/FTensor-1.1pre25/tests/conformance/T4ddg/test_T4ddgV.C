#include <iostream>
#include "../../../FTensor.h"
#include "../test_for_zero.h"
using namespace FTensor;
using namespace std;

void test_T4ddgV(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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

  /* T4_ddg(1,1,i,j)=T2s(i,j) */

  Tensor1<double,3> t1_3;
  t1_3(i)=t1_1(j)*t2_1(i,j);

  t4ddg_1(N0,N0,i,j)=(t1_1(i)^t1_1(j));
  t4ddg_1(N0,N1,i,j)=t2s_1(i,j);
  t4ddg_1(N0,N2,i,j)=(t1_3(i)^t1_3(j));
  t4ddg_1(N1,N1,i,j)=t2s_2(i,j);
  t4ddg_1(N1,N2,i,j)=(t1_2(i)^t1_2(j));
  t4ddg_1(N2,N2,i,j)=t2s_3(i,j);

  test_for_zero(t4ddg_1(0,0,0,0) - (t1_1(0)*t1_1(0))
		,"T4ddg(N,N,i,j)=T2s(i,j)(0,0,0,0)");
  test_for_zero(t4ddg_1(0,0,0,1) - (t1_1(0)*t1_1(1))
		,"T4ddg(N,N,i,j)=T2s(i,j)(0,0,0,1)");
  test_for_zero(t4ddg_1(0,0,0,2) - (t1_1(0)*t1_1(2))
		,"T4ddg(N,N,i,j)=T2s(i,j)(0,0,0,2)");
  test_for_zero(t4ddg_1(0,0,1,0) - (t1_1(1)*t1_1(0))
		,"T4ddg(N,N,i,j)=T2s(i,j)(0,0,1,0)");
  test_for_zero(t4ddg_1(0,0,1,1) - (t1_1(1)*t1_1(1))
		,"T4ddg(N,N,i,j)=T2s(i,j)(0,0,1,1)");
  test_for_zero(t4ddg_1(0,0,1,2) - (t1_1(1)*t1_1(2))
		,"T4ddg(N,N,i,j)=T2s(i,j)(0,0,1,2)");
  test_for_zero(t4ddg_1(0,0,2,0) - (t1_1(2)*t1_1(0))
		,"T4ddg(N,N,i,j)=T2s(i,j)(0,0,2,0)");
  test_for_zero(t4ddg_1(0,0,2,1) - (t1_1(2)*t1_1(1))
		,"T4ddg(N,N,i,j)=T2s(i,j)(0,0,2,1)");
  test_for_zero(t4ddg_1(0,0,2,2) - (t1_1(2)*t1_1(2))
		,"T4ddg(N,N,i,j)=T2s(i,j)(0,0,2,2)");
  test_for_zero(t4ddg_1(0,1,0,0) - (t2s_1(0,0))
		,"T4ddg(N,N,i,j)=T2s(i,j)(0,1,0,0)");
  test_for_zero(t4ddg_1(0,1,0,1) - (t2s_1(0,1))
		,"T4ddg(N,N,i,j)=T2s(i,j)(0,1,0,1)");
  test_for_zero(t4ddg_1(0,1,0,2) - (t2s_1(0,2))
		,"T4ddg(N,N,i,j)=T2s(i,j)(0,1,0,2)");
  test_for_zero(t4ddg_1(0,1,1,0) - (t2s_1(1,0))
		,"T4ddg(N,N,i,j)=T2s(i,j)(0,1,1,0)");
  test_for_zero(t4ddg_1(0,1,1,1) - (t2s_1(1,1))
		,"T4ddg(N,N,i,j)=T2s(i,j)(0,1,1,1)");
  test_for_zero(t4ddg_1(0,1,1,2) - (t2s_1(1,2))
		,"T4ddg(N,N,i,j)=T2s(i,j)(0,1,1,2)");
  test_for_zero(t4ddg_1(0,1,2,0) - (t2s_1(2,0))
		,"T4ddg(N,N,i,j)=T2s(i,j)(0,1,2,0)");
  test_for_zero(t4ddg_1(0,1,2,1) - (t2s_1(2,1))
		,"T4ddg(N,N,i,j)=T2s(i,j)(0,1,2,1)");
  test_for_zero(t4ddg_1(0,1,2,2) - (t2s_1(2,2))
		,"T4ddg(N,N,i,j)=T2s(i,j)(0,1,2,2)");
  test_for_zero((t4ddg_1(0,2,0,0) - (t1_3(0)*t1_3(0)))/t4ddg_1(0,2,0,0)
		,"T4ddg(N,N,i,j)=T2s(i,j)(0,2,0,0)");
  test_for_zero((t4ddg_1(0,2,0,1) - (t1_3(0)*t1_3(1)))/t4ddg_1(0,2,0,1)
		,"T4ddg(N,N,i,j)=T2s(i,j)(0,2,0,1)");
  test_for_zero((t4ddg_1(0,2,0,2) - (t1_3(0)*t1_3(2)))/t4ddg_1(0,2,0,2)
		,"T4ddg(N,N,i,j)=T2s(i,j)(0,2,0,2)");
  test_for_zero((t4ddg_1(0,2,1,0) - (t1_3(1)*t1_3(0)))/t4ddg_1(0,2,1,0)
		,"T4ddg(N,N,i,j)=T2s(i,j)(0,2,1,0)");
  test_for_zero((t4ddg_1(0,2,1,1) - (t1_3(1)*t1_3(1)))/t4ddg_1(0,2,1,1)
		,"T4ddg(N,N,i,j)=T2s(i,j)(0,2,1,1)");
  test_for_zero((t4ddg_1(0,2,1,2) - (t1_3(1)*t1_3(2)))/t4ddg_1(0,2,1,2)
		,"T4ddg(N,N,i,j)=T2s(i,j)(0,2,1,2)");
  test_for_zero((t4ddg_1(0,2,2,0) - (t1_3(2)*t1_3(0)))/t4ddg_1(0,2,2,0)
		,"T4ddg(N,N,i,j)=T2s(i,j)(0,2,2,0)");
  test_for_zero((t4ddg_1(0,2,2,1) - (t1_3(2)*t1_3(1)))/t4ddg_1(0,2,2,1)
		,"T4ddg(N,N,i,j)=T2s(i,j)(0,2,2,1)");
  test_for_zero((t4ddg_1(0,2,2,2) - (t1_3(2)*t1_3(2)))/t4ddg_1(0,2,2,2)
		,"T4ddg(N,N,i,j)=T2s(i,j)(0,2,2,2)");
  test_for_zero(t4ddg_1(1,0,0,0) - (t2s_1(0,0))
		,"T4ddg(N,N,i,j)=T2s(i,j)(1,0,0,0)");
  test_for_zero(t4ddg_1(1,0,0,1) - (t2s_1(0,1))
		,"T4ddg(N,N,i,j)=T2s(i,j)(1,0,0,1)");
  test_for_zero(t4ddg_1(1,0,0,2) - (t2s_1(0,2))
		,"T4ddg(N,N,i,j)=T2s(i,j)(1,0,0,2)");
  test_for_zero(t4ddg_1(1,0,1,0) - (t2s_1(1,0))
		,"T4ddg(N,N,i,j)=T2s(i,j)(1,0,1,0)");
  test_for_zero(t4ddg_1(1,0,1,1) - (t2s_1(1,1))
		,"T4ddg(N,N,i,j)=T2s(i,j)(1,0,1,1)");
  test_for_zero(t4ddg_1(1,0,1,2) - (t2s_1(1,2))
		,"T4ddg(N,N,i,j)=T2s(i,j)(1,0,1,2)");
  test_for_zero(t4ddg_1(1,0,2,0) - (t2s_1(2,0))
		,"T4ddg(N,N,i,j)=T2s(i,j)(1,0,2,0)");
  test_for_zero(t4ddg_1(1,0,2,1) - (t2s_1(2,1))
		,"T4ddg(N,N,i,j)=T2s(i,j)(1,0,2,1)");
  test_for_zero(t4ddg_1(1,0,2,2) - (t2s_1(2,2))
		,"T4ddg(N,N,i,j)=T2s(i,j)(1,0,2,2)");
  test_for_zero(t4ddg_1(1,1,0,0) - (t2s_2(0,0))
		,"T4ddg(N,N,i,j)=T2s(i,j)(1,1,0,0)");
  test_for_zero(t4ddg_1(1,1,0,1) - (t2s_2(0,1))
		,"T4ddg(N,N,i,j)=T2s(i,j)(1,1,0,1)");
  test_for_zero(t4ddg_1(1,1,0,2) - (t2s_2(0,2))
		,"T4ddg(N,N,i,j)=T2s(i,j)(1,1,0,2)");
  test_for_zero(t4ddg_1(1,1,1,0) - (t2s_2(1,0))
		,"T4ddg(N,N,i,j)=T2s(i,j)(1,1,1,0)");
  test_for_zero(t4ddg_1(1,1,1,1) - (t2s_2(1,1))
		,"T4ddg(N,N,i,j)=T2s(i,j)(1,1,1,1)");
  test_for_zero(t4ddg_1(1,1,1,2) - (t2s_2(1,2))
		,"T4ddg(N,N,i,j)=T2s(i,j)(1,1,1,2)");
  test_for_zero(t4ddg_1(1,1,2,0) - (t2s_2(2,0))
		,"T4ddg(N,N,i,j)=T2s(i,j)(1,1,2,0)");
  test_for_zero(t4ddg_1(1,1,2,1) - (t2s_2(2,1))
		,"T4ddg(N,N,i,j)=T2s(i,j)(1,1,2,1)");
  test_for_zero(t4ddg_1(1,1,2,2) - (t2s_2(2,2))
		,"T4ddg(N,N,i,j)=T2s(i,j)(1,1,2,2)");
  test_for_zero(t4ddg_1(1,2,0,0) - (t1_2(0)*t1_2(0))
		,"T4ddg(N,N,i,j)=T2s(i,j)(1,2,0,0)");
  test_for_zero(t4ddg_1(1,2,0,1) - (t1_2(0)*t1_2(1))
		,"T4ddg(N,N,i,j)=T2s(i,j)(1,2,0,1)");
  test_for_zero(t4ddg_1(1,2,0,2) - (t1_2(0)*t1_2(2))
		,"T4ddg(N,N,i,j)=T2s(i,j)(1,2,0,2)");
  test_for_zero(t4ddg_1(1,2,1,0) - (t1_2(1)*t1_2(0))
		,"T4ddg(N,N,i,j)=T2s(i,j)(1,2,1,0)");
  test_for_zero(t4ddg_1(1,2,1,1) - (t1_2(1)*t1_2(1))
		,"T4ddg(N,N,i,j)=T2s(i,j)(1,2,1,1)");
  test_for_zero(t4ddg_1(1,2,1,2) - (t1_2(1)*t1_2(2))
		,"T4ddg(N,N,i,j)=T2s(i,j)(1,2,1,2)");
  test_for_zero(t4ddg_1(1,2,2,0) - (t1_2(2)*t1_2(0))
		,"T4ddg(N,N,i,j)=T2s(i,j)(1,2,2,0)");
  test_for_zero(t4ddg_1(1,2,2,1) - (t1_2(2)*t1_2(1))
		,"T4ddg(N,N,i,j)=T2s(i,j)(1,2,2,1)");
  test_for_zero(t4ddg_1(1,2,2,2) - (t1_2(2)*t1_2(2))
		,"T4ddg(N,N,i,j)=T2s(i,j)(1,2,2,2)");
  test_for_zero((t4ddg_1(2,0,0,0) - (t1_3(0)*t1_3(0)))/t4ddg_1(2,0,0,0)
		,"T4ddg(N,N,i,j)=T2s(i,j)(2,0,0,0)");
  test_for_zero((t4ddg_1(2,0,0,1) - (t1_3(0)*t1_3(1)))/t4ddg_1(2,0,0,1)
		,"T4ddg(N,N,i,j)=T2s(i,j)(2,0,0,1)");
  test_for_zero((t4ddg_1(2,0,0,2) - (t1_3(0)*t1_3(2)))/t4ddg_1(2,0,0,2)
		,"T4ddg(N,N,i,j)=T2s(i,j)(2,0,0,2)");
  test_for_zero((t4ddg_1(2,0,1,0) - (t1_3(1)*t1_3(0)))/t4ddg_1(2,0,1,0)
		,"T4ddg(N,N,i,j)=T2s(i,j)(2,0,1,0)");
  test_for_zero((t4ddg_1(2,0,1,1) - (t1_3(1)*t1_3(1)))/t4ddg_1(2,0,1,1)
		,"T4ddg(N,N,i,j)=T2s(i,j)(2,0,1,1)");
  test_for_zero((t4ddg_1(2,0,1,2) - (t1_3(1)*t1_3(2)))/t4ddg_1(2,0,1,2)
		,"T4ddg(N,N,i,j)=T2s(i,j)(2,0,1,2)");
  test_for_zero((t4ddg_1(2,0,2,0) - (t1_3(2)*t1_3(0)))/t4ddg_1(2,0,2,0)
		,"T4ddg(N,N,i,j)=T2s(i,j)(2,0,2,0)");
  test_for_zero((t4ddg_1(2,0,2,1) - (t1_3(2)*t1_3(1)))/t4ddg_1(2,0,2,1)
		,"T4ddg(N,N,i,j)=T2s(i,j)(2,0,2,1)");
  test_for_zero((t4ddg_1(2,0,2,2) - (t1_3(2)*t1_3(2)))/t4ddg_1(2,0,2,2)
		,"T4ddg(N,N,i,j)=T2s(i,j)(2,0,2,2)");
  test_for_zero(t4ddg_1(2,1,0,0) - (t1_2(0)*t1_2(0))
		,"T4ddg(N,N,i,j)=T2s(i,j)(2,1,0,0)");
  test_for_zero(t4ddg_1(2,1,0,1) - (t1_2(0)*t1_2(1))
		,"T4ddg(N,N,i,j)=T2s(i,j)(2,1,0,1)");
  test_for_zero(t4ddg_1(2,1,0,2) - (t1_2(0)*t1_2(2))
		,"T4ddg(N,N,i,j)=T2s(i,j)(2,1,0,2)");
  test_for_zero(t4ddg_1(2,1,1,0) - (t1_2(1)*t1_2(0))
		,"T4ddg(N,N,i,j)=T2s(i,j)(2,1,1,0)");
  test_for_zero(t4ddg_1(2,1,1,1) - (t1_2(1)*t1_2(1))
		,"T4ddg(N,N,i,j)=T2s(i,j)(2,1,1,1)");
  test_for_zero(t4ddg_1(2,1,1,2) - (t1_2(1)*t1_2(2))
		,"T4ddg(N,N,i,j)=T2s(i,j)(2,1,1,2)");
  test_for_zero(t4ddg_1(2,1,2,0) - (t1_2(2)*t1_2(0))
		,"T4ddg(N,N,i,j)=T2s(i,j)(2,1,2,0)");
  test_for_zero(t4ddg_1(2,1,2,1) - (t1_2(2)*t1_2(1))
		,"T4ddg(N,N,i,j)=T2s(i,j)(2,1,2,1)");
  test_for_zero(t4ddg_1(2,1,2,2) - (t1_2(2)*t1_2(2))
		,"T4ddg(N,N,i,j)=T2s(i,j)(2,1,2,2)");
  test_for_zero(t4ddg_1(2,2,0,0) - (t2s_3(0,0))
		,"T4ddg(N,N,i,j)=T2s(i,j)(2,2,0,0)");
  test_for_zero(t4ddg_1(2,2,0,1) - (t2s_3(0,1))
		,"T4ddg(N,N,i,j)=T2s(i,j)(2,2,0,1)");
  test_for_zero(t4ddg_1(2,2,0,2) - (t2s_3(0,2))
		,"T4ddg(N,N,i,j)=T2s(i,j)(2,2,0,2)");
  test_for_zero(t4ddg_1(2,2,1,0) - (t2s_3(1,0))
		,"T4ddg(N,N,i,j)=T2s(i,j)(2,2,1,0)");
  test_for_zero(t4ddg_1(2,2,1,1) - (t2s_3(1,1))
		,"T4ddg(N,N,i,j)=T2s(i,j)(2,2,1,1)");
  test_for_zero(t4ddg_1(2,2,1,2) - (t2s_3(1,2))
		,"T4ddg(N,N,i,j)=T2s(i,j)(2,2,1,2)");
  test_for_zero(t4ddg_1(2,2,2,0) - (t2s_3(2,0))
		,"T4ddg(N,N,i,j)=T2s(i,j)(2,2,2,0)");
  test_for_zero(t4ddg_1(2,2,2,1) - (t2s_3(2,1))
		,"T4ddg(N,N,i,j)=T2s(i,j)(2,2,2,1)");
  test_for_zero(t4ddg_1(2,2,2,2) - (t2s_3(2,2))
		,"T4ddg(N,N,i,j)=T2s(i,j)(2,2,2,2)");

  t4ddg_2(i,j,k,l)=13*t4ddg_1(i,j,k,l);
  test_for_zero(t4ddg_2(0,0,0,0) - 13*(t1_1(0)*t1_1(0))
		,"T*T4ddg(0,0,0,0)");
  test_for_zero(t4ddg_2(0,0,0,1) - 13*(t1_1(0)*t1_1(1))
		,"T*T4ddg(0,0,0,1)");
  test_for_zero(t4ddg_2(0,0,0,2) - 13*(t1_1(0)*t1_1(2))
		,"T*T4ddg(0,0,0,2)");
  test_for_zero(t4ddg_2(0,0,1,0) - 13*(t1_1(1)*t1_1(0))
		,"T*T4ddg(0,0,1,0)");
  test_for_zero(t4ddg_2(0,0,1,1) - 13*(t1_1(1)*t1_1(1))
		,"T*T4ddg(0,0,1,1)");
  test_for_zero(t4ddg_2(0,0,1,2) - 13*(t1_1(1)*t1_1(2))
		,"T*T4ddg(0,0,1,2)");
  test_for_zero(t4ddg_2(0,0,2,0) - 13*(t1_1(2)*t1_1(0))
		,"T*T4ddg(0,0,2,0)");
  test_for_zero(t4ddg_2(0,0,2,1) - 13*(t1_1(2)*t1_1(1))
		,"T*T4ddg(0,0,2,1)");
  test_for_zero(t4ddg_2(0,0,2,2) - 13*(t1_1(2)*t1_1(2))
		,"T*T4ddg(0,0,2,2)");
  test_for_zero(t4ddg_2(0,1,0,0) - 13*(t2s_1(0,0))
		,"T*T4ddg(0,1,0,0)");
  test_for_zero(t4ddg_2(0,1,0,1) - 13*(t2s_1(0,1))
		,"T*T4ddg(0,1,0,1)");
  test_for_zero(t4ddg_2(0,1,0,2) - 13*(t2s_1(0,2))
		,"T*T4ddg(0,1,0,2)");
  test_for_zero(t4ddg_2(0,1,1,0) - 13*(t2s_1(1,0))
		,"T*T4ddg(0,1,1,0)");
  test_for_zero(t4ddg_2(0,1,1,1) - 13*(t2s_1(1,1))
		,"T*T4ddg(0,1,1,1)");
  test_for_zero(t4ddg_2(0,1,1,2) - 13*(t2s_1(1,2))
		,"T*T4ddg(0,1,1,2)");
  test_for_zero(t4ddg_2(0,1,2,0) - 13*(t2s_1(2,0))
		,"T*T4ddg(0,1,2,0)");
  test_for_zero(t4ddg_2(0,1,2,1) - 13*(t2s_1(2,1))
		,"T*T4ddg(0,1,2,1)");
  test_for_zero(t4ddg_2(0,1,2,2) - 13*(t2s_1(2,2))
		,"T*T4ddg(0,1,2,2)");
  test_for_zero((t4ddg_2(0,2,0,0) - 13*(t1_3(0)*t1_3(0)))/t4ddg_2(0,2,0,0)
		,"T*T4ddg(0,2,0,0)");
  test_for_zero((t4ddg_2(0,2,0,1) - 13*(t1_3(0)*t1_3(1)))/t4ddg_2(0,2,0,1)
		,"T*T4ddg(0,2,0,1)");
  test_for_zero((t4ddg_2(0,2,0,2) - 13*(t1_3(0)*t1_3(2)))/t4ddg_2(0,2,0,2)
		,"T*T4ddg(0,2,0,2)");
  test_for_zero((t4ddg_2(0,2,1,0) - 13*(t1_3(1)*t1_3(0)))/t4ddg_2(0,2,1,0)
		,"T*T4ddg(0,2,1,0)");
  test_for_zero((t4ddg_2(0,2,1,1) - 13*(t1_3(1)*t1_3(1)))/t4ddg_2(0,2,1,1)
		,"T*T4ddg(0,2,1,1)");
  test_for_zero((t4ddg_2(0,2,1,2) - 13*(t1_3(1)*t1_3(2)))/t4ddg_2(0,2,1,2)
		,"T*T4ddg(0,2,1,2)");
  test_for_zero((t4ddg_2(0,2,2,0) - 13*(t1_3(2)*t1_3(0)))/t4ddg_2(0,2,2,0)
		,"T*T4ddg(0,2,2,0)");
  test_for_zero((t4ddg_2(0,2,2,1) - 13*(t1_3(2)*t1_3(1)))/t4ddg_2(0,2,2,1)
		,"T*T4ddg(0,2,2,1)");
  test_for_zero((t4ddg_2(0,2,2,2) - 13*(t1_3(2)*t1_3(2)))/t4ddg_2(0,2,2,2)
		,"T*T4ddg(0,2,2,2)");
  test_for_zero(t4ddg_2(1,0,0,0) - 13*(t2s_1(0,0))
		,"T*T4ddg(1,0,0,0)");
  test_for_zero(t4ddg_2(1,0,0,1) - 13*(t2s_1(0,1))
		,"T*T4ddg(1,0,0,1)");
  test_for_zero(t4ddg_2(1,0,0,2) - 13*(t2s_1(0,2))
		,"T*T4ddg(1,0,0,2)");
  test_for_zero(t4ddg_2(1,0,1,0) - 13*(t2s_1(1,0))
		,"T*T4ddg(1,0,1,0)");
  test_for_zero(t4ddg_2(1,0,1,1) - 13*(t2s_1(1,1))
		,"T*T4ddg(1,0,1,1)");
  test_for_zero(t4ddg_2(1,0,1,2) - 13*(t2s_1(1,2))
		,"T*T4ddg(1,0,1,2)");
  test_for_zero(t4ddg_2(1,0,2,0) - 13*(t2s_1(2,0))
		,"T*T4ddg(1,0,2,0)");
  test_for_zero(t4ddg_2(1,0,2,1) - 13*(t2s_1(2,1))
		,"T*T4ddg(1,0,2,1)");
  test_for_zero(t4ddg_2(1,0,2,2) - 13*(t2s_1(2,2))
		,"T*T4ddg(1,0,2,2)");
  test_for_zero(t4ddg_2(1,1,0,0) - 13*(t2s_2(0,0))
		,"T*T4ddg(1,1,0,0)");
  test_for_zero(t4ddg_2(1,1,0,1) - 13*(t2s_2(0,1))
		,"T*T4ddg(1,1,0,1)");
  test_for_zero(t4ddg_2(1,1,0,2) - 13*(t2s_2(0,2))
		,"T*T4ddg(1,1,0,2)");
  test_for_zero(t4ddg_2(1,1,1,0) - 13*(t2s_2(1,0))
		,"T*T4ddg(1,1,1,0)");
  test_for_zero(t4ddg_2(1,1,1,1) - 13*(t2s_2(1,1))
		,"T*T4ddg(1,1,1,1)");
  test_for_zero(t4ddg_2(1,1,1,2) - 13*(t2s_2(1,2))
		,"T*T4ddg(1,1,1,2)");
  test_for_zero(t4ddg_2(1,1,2,0) - 13*(t2s_2(2,0))
		,"T*T4ddg(1,1,2,0)");
  test_for_zero(t4ddg_2(1,1,2,1) - 13*(t2s_2(2,1))
		,"T*T4ddg(1,1,2,1)");
  test_for_zero(t4ddg_2(1,1,2,2) - 13*(t2s_2(2,2))
		,"T*T4ddg(1,1,2,2)");
  test_for_zero(t4ddg_2(1,2,0,0) - 13*(t1_2(0)*t1_2(0))
		,"T*T4ddg(1,2,0,0)");
  test_for_zero(t4ddg_2(1,2,0,1) - 13*(t1_2(0)*t1_2(1))
		,"T*T4ddg(1,2,0,1)");
  test_for_zero(t4ddg_2(1,2,0,2) - 13*(t1_2(0)*t1_2(2))
		,"T*T4ddg(1,2,0,2)");
  test_for_zero(t4ddg_2(1,2,1,0) - 13*(t1_2(1)*t1_2(0))
		,"T*T4ddg(1,2,1,0)");
  test_for_zero(t4ddg_2(1,2,1,1) - 13*(t1_2(1)*t1_2(1))
		,"T*T4ddg(1,2,1,1)");
  test_for_zero(t4ddg_2(1,2,1,2) - 13*(t1_2(1)*t1_2(2))
		,"T*T4ddg(1,2,1,2)");
  test_for_zero(t4ddg_2(1,2,2,0) - 13*(t1_2(2)*t1_2(0))
		,"T*T4ddg(1,2,2,0)");
  test_for_zero(t4ddg_2(1,2,2,1) - 13*(t1_2(2)*t1_2(1))
		,"T*T4ddg(1,2,2,1)");
  test_for_zero(t4ddg_2(1,2,2,2) - 13*(t1_2(2)*t1_2(2))
		,"T*T4ddg(1,2,2,2)");
  test_for_zero((t4ddg_2(2,0,0,0) - 13*(t1_3(0)*t1_3(0)))/t4ddg_2(2,0,0,0)
		,"T*T4ddg(2,0,0,0)");
  test_for_zero((t4ddg_2(2,0,0,1) - 13*(t1_3(0)*t1_3(1)))/t4ddg_2(2,0,0,1)
		,"T*T4ddg(2,0,0,1)");
  test_for_zero((t4ddg_2(2,0,0,2) - 13*(t1_3(0)*t1_3(2)))/t4ddg_2(2,0,0,2)
		,"T*T4ddg(2,0,0,2)");
  test_for_zero((t4ddg_2(2,0,1,0) - 13*(t1_3(1)*t1_3(0)))/t4ddg_2(2,0,1,0)
		,"T*T4ddg(2,0,1,0)");
  test_for_zero((t4ddg_2(2,0,1,1) - 13*(t1_3(1)*t1_3(1)))/t4ddg_2(2,0,1,1)
		,"T*T4ddg(2,0,1,1)");
  test_for_zero((t4ddg_2(2,0,1,2) - 13*(t1_3(1)*t1_3(2)))/t4ddg_2(2,0,1,2)
		,"T*T4ddg(2,0,1,2)");
  test_for_zero((t4ddg_2(2,0,2,0) - 13*(t1_3(2)*t1_3(0)))/t4ddg_2(2,0,2,0)
		,"T*T4ddg(2,0,2,0)");
  test_for_zero((t4ddg_2(2,0,2,1) - 13*(t1_3(2)*t1_3(1)))/t4ddg_2(2,0,2,1)
		,"T*T4ddg(2,0,2,1)");
  test_for_zero((t4ddg_2(2,0,2,2) - 13*(t1_3(2)*t1_3(2)))/t4ddg_2(2,0,2,2)
		,"T*T4ddg(2,0,2,2)");
  test_for_zero(t4ddg_2(2,1,0,0) - 13*(t1_2(0)*t1_2(0))
		,"T*T4ddg(2,1,0,0)");
  test_for_zero(t4ddg_2(2,1,0,1) - 13*(t1_2(0)*t1_2(1))
		,"T*T4ddg(2,1,0,1)");
  test_for_zero(t4ddg_2(2,1,0,2) - 13*(t1_2(0)*t1_2(2))
		,"T*T4ddg(2,1,0,2)");
  test_for_zero(t4ddg_2(2,1,1,0) - 13*(t1_2(1)*t1_2(0))
		,"T*T4ddg(2,1,1,0)");
  test_for_zero(t4ddg_2(2,1,1,1) - 13*(t1_2(1)*t1_2(1))
		,"T*T4ddg(2,1,1,1)");
  test_for_zero(t4ddg_2(2,1,1,2) - 13*(t1_2(1)*t1_2(2))
		,"T*T4ddg(2,1,1,2)");
  test_for_zero(t4ddg_2(2,1,2,0) - 13*(t1_2(2)*t1_2(0))
		,"T*T4ddg(2,1,2,0)");
  test_for_zero(t4ddg_2(2,1,2,1) - 13*(t1_2(2)*t1_2(1))
		,"T*T4ddg(2,1,2,1)");
  test_for_zero(t4ddg_2(2,1,2,2) - 13*(t1_2(2)*t1_2(2))
		,"T*T4ddg(2,1,2,2)");
  test_for_zero(t4ddg_2(2,2,0,0) - 13*(t2s_3(0,0))
		,"T*T4ddg(2,2,0,0)");
  test_for_zero(t4ddg_2(2,2,0,1) - 13*(t2s_3(0,1))
		,"T*T4ddg(2,2,0,1)");
  test_for_zero(t4ddg_2(2,2,0,2) - 13*(t2s_3(0,2))
		,"T*T4ddg(2,2,0,2)");
  test_for_zero(t4ddg_2(2,2,1,0) - 13*(t2s_3(1,0))
		,"T*T4ddg(2,2,1,0)");
  test_for_zero(t4ddg_2(2,2,1,1) - 13*(t2s_3(1,1))
		,"T*T4ddg(2,2,1,1)");
  test_for_zero(t4ddg_2(2,2,1,2) - 13*(t2s_3(1,2))
		,"T*T4ddg(2,2,1,2)");
  test_for_zero(t4ddg_2(2,2,2,0) - 13*(t2s_3(2,0))
		,"T*T4ddg(2,2,2,0)");
  test_for_zero(t4ddg_2(2,2,2,1) - 13*(t2s_3(2,1))
		,"T*T4ddg(2,2,2,1)");
  test_for_zero(t4ddg_2(2,2,2,2) - 13*(t2s_3(2,2))
		,"T*T4ddg(2,2,2,2)");


  t4ddg_2(i,j,k,l)=t4ddg_1(i,j,k,l)*7;
  test_for_zero(t4ddg_2(0,0,0,0) - 7*(t1_1(0)*t1_1(0))
		,"T4ddg*T(0,0,0,0)");
  test_for_zero(t4ddg_2(0,0,0,1) - 7*(t1_1(0)*t1_1(1))
		,"T4ddg*T(0,0,0,1)");
  test_for_zero(t4ddg_2(0,0,0,2) - 7*(t1_1(0)*t1_1(2))
		,"T4ddg*T(0,0,0,2)");
  test_for_zero(t4ddg_2(0,0,1,0) - 7*(t1_1(1)*t1_1(0))
		,"T4ddg*T(0,0,1,0)");
  test_for_zero(t4ddg_2(0,0,1,1) - 7*(t1_1(1)*t1_1(1))
		,"T4ddg*T(0,0,1,1)");
  test_for_zero(t4ddg_2(0,0,1,2) - 7*(t1_1(1)*t1_1(2))
		,"T4ddg*T(0,0,1,2)");
  test_for_zero(t4ddg_2(0,0,2,0) - 7*(t1_1(2)*t1_1(0))
		,"T4ddg*T(0,0,2,0)");
  test_for_zero(t4ddg_2(0,0,2,1) - 7*(t1_1(2)*t1_1(1))
		,"T4ddg*T(0,0,2,1)");
  test_for_zero(t4ddg_2(0,0,2,2) - 7*(t1_1(2)*t1_1(2))
		,"T4ddg*T(0,0,2,2)");
  test_for_zero(t4ddg_2(0,1,0,0) - 7*(t2s_1(0,0))
		,"T4ddg*T(0,1,0,0)");
  test_for_zero(t4ddg_2(0,1,0,1) - 7*(t2s_1(0,1))
		,"T4ddg*T(0,1,0,1)");
  test_for_zero(t4ddg_2(0,1,0,2) - 7*(t2s_1(0,2))
		,"T4ddg*T(0,1,0,2)");
  test_for_zero(t4ddg_2(0,1,1,0) - 7*(t2s_1(1,0))
		,"T4ddg*T(0,1,1,0)");
  test_for_zero(t4ddg_2(0,1,1,1) - 7*(t2s_1(1,1))
		,"T4ddg*T(0,1,1,1)");
  test_for_zero(t4ddg_2(0,1,1,2) - 7*(t2s_1(1,2))
		,"T4ddg*T(0,1,1,2)");
  test_for_zero(t4ddg_2(0,1,2,0) - 7*(t2s_1(2,0))
		,"T4ddg*T(0,1,2,0)");
  test_for_zero(t4ddg_2(0,1,2,1) - 7*(t2s_1(2,1))
		,"T4ddg*T(0,1,2,1)");
  test_for_zero(t4ddg_2(0,1,2,2) - 7*(t2s_1(2,2))
		,"T4ddg*T(0,1,2,2)");
  test_for_zero((t4ddg_2(0,2,0,0) - 7*(t1_3(0)*t1_3(0)))/t4ddg_2(0,2,0,0)
		,"T4ddg*T(0,2,0,0)");
  test_for_zero((t4ddg_2(0,2,0,1) - 7*(t1_3(0)*t1_3(1)))/t4ddg_2(0,2,0,1)
		,"T4ddg*T(0,2,0,1)");
  test_for_zero((t4ddg_2(0,2,0,2) - 7*(t1_3(0)*t1_3(2)))/t4ddg_2(0,2,0,2)
		,"T4ddg*T(0,2,0,2)");
  test_for_zero((t4ddg_2(0,2,1,0) - 7*(t1_3(1)*t1_3(0)))/t4ddg_2(0,2,1,0)
		,"T4ddg*T(0,2,1,0)");
  test_for_zero((t4ddg_2(0,2,1,1) - 7*(t1_3(1)*t1_3(1)))/t4ddg_2(0,2,1,1)
		,"T4ddg*T(0,2,1,1)");
  test_for_zero((t4ddg_2(0,2,1,2) - 7*(t1_3(1)*t1_3(2)))/t4ddg_2(0,2,1,2)
		,"T4ddg*T(0,2,1,2)");
  test_for_zero((t4ddg_2(0,2,2,0) - 7*(t1_3(2)*t1_3(0)))/t4ddg_2(0,2,2,0)
		,"T4ddg*T(0,2,2,0)");
  test_for_zero((t4ddg_2(0,2,2,1) - 7*(t1_3(2)*t1_3(1)))/t4ddg_2(0,2,2,1)
		,"T4ddg*T(0,2,2,1)");
  test_for_zero((t4ddg_2(0,2,2,2) - 7*(t1_3(2)*t1_3(2)))/t4ddg_2(0,2,2,2)
		,"T4ddg*T(0,2,2,2)");
  test_for_zero(t4ddg_2(1,0,0,0) - 7*(t2s_1(0,0))
		,"T4ddg*T(1,0,0,0)");
  test_for_zero(t4ddg_2(1,0,0,1) - 7*(t2s_1(0,1))
		,"T4ddg*T(1,0,0,1)");
  test_for_zero(t4ddg_2(1,0,0,2) - 7*(t2s_1(0,2))
		,"T4ddg*T(1,0,0,2)");
  test_for_zero(t4ddg_2(1,0,1,0) - 7*(t2s_1(1,0))
		,"T4ddg*T(1,0,1,0)");
  test_for_zero(t4ddg_2(1,0,1,1) - 7*(t2s_1(1,1))
		,"T4ddg*T(1,0,1,1)");
  test_for_zero(t4ddg_2(1,0,1,2) - 7*(t2s_1(1,2))
		,"T4ddg*T(1,0,1,2)");
  test_for_zero(t4ddg_2(1,0,2,0) - 7*(t2s_1(2,0))
		,"T4ddg*T(1,0,2,0)");
  test_for_zero(t4ddg_2(1,0,2,1) - 7*(t2s_1(2,1))
		,"T4ddg*T(1,0,2,1)");
  test_for_zero(t4ddg_2(1,0,2,2) - 7*(t2s_1(2,2))
		,"T4ddg*T(1,0,2,2)");
  test_for_zero(t4ddg_2(1,1,0,0) - 7*(t2s_2(0,0))
		,"T4ddg*T(1,1,0,0)");
  test_for_zero(t4ddg_2(1,1,0,1) - 7*(t2s_2(0,1))
		,"T4ddg*T(1,1,0,1)");
  test_for_zero(t4ddg_2(1,1,0,2) - 7*(t2s_2(0,2))
		,"T4ddg*T(1,1,0,2)");
  test_for_zero(t4ddg_2(1,1,1,0) - 7*(t2s_2(1,0))
		,"T4ddg*T(1,1,1,0)");
  test_for_zero(t4ddg_2(1,1,1,1) - 7*(t2s_2(1,1))
		,"T4ddg*T(1,1,1,1)");
  test_for_zero(t4ddg_2(1,1,1,2) - 7*(t2s_2(1,2))
		,"T4ddg*T(1,1,1,2)");
  test_for_zero(t4ddg_2(1,1,2,0) - 7*(t2s_2(2,0))
		,"T4ddg*T(1,1,2,0)");
  test_for_zero(t4ddg_2(1,1,2,1) - 7*(t2s_2(2,1))
		,"T4ddg*T(1,1,2,1)");
  test_for_zero(t4ddg_2(1,1,2,2) - 7*(t2s_2(2,2))
		,"T4ddg*T(1,1,2,2)");
  test_for_zero(t4ddg_2(1,2,0,0) - 7*(t1_2(0)*t1_2(0))
		,"T4ddg*T(1,2,0,0)");
  test_for_zero(t4ddg_2(1,2,0,1) - 7*(t1_2(0)*t1_2(1))
		,"T4ddg*T(1,2,0,1)");
  test_for_zero(t4ddg_2(1,2,0,2) - 7*(t1_2(0)*t1_2(2))
		,"T4ddg*T(1,2,0,2)");
  test_for_zero(t4ddg_2(1,2,1,0) - 7*(t1_2(1)*t1_2(0))
		,"T4ddg*T(1,2,1,0)");
  test_for_zero(t4ddg_2(1,2,1,1) - 7*(t1_2(1)*t1_2(1))
		,"T4ddg*T(1,2,1,1)");
  test_for_zero(t4ddg_2(1,2,1,2) - 7*(t1_2(1)*t1_2(2))
		,"T4ddg*T(1,2,1,2)");
  test_for_zero(t4ddg_2(1,2,2,0) - 7*(t1_2(2)*t1_2(0))
		,"T4ddg*T(1,2,2,0)");
  test_for_zero(t4ddg_2(1,2,2,1) - 7*(t1_2(2)*t1_2(1))
		,"T4ddg*T(1,2,2,1)");
  test_for_zero(t4ddg_2(1,2,2,2) - 7*(t1_2(2)*t1_2(2))
		,"T4ddg*T(1,2,2,2)");
  test_for_zero((t4ddg_2(2,0,0,0) - 7*(t1_3(0)*t1_3(0)))/t4ddg_2(2,0,0,0)
		,"T4ddg*T(2,0,0,0)");
  test_for_zero((t4ddg_2(2,0,0,1) - 7*(t1_3(0)*t1_3(1)))/t4ddg_2(2,0,0,1)
		,"T4ddg*T(2,0,0,1)");
  test_for_zero((t4ddg_2(2,0,0,2) - 7*(t1_3(0)*t1_3(2)))/t4ddg_2(2,0,0,2)
		,"T4ddg*T(2,0,0,2)");
  test_for_zero((t4ddg_2(2,0,1,0) - 7*(t1_3(1)*t1_3(0)))/t4ddg_2(2,0,1,0)
		,"T4ddg*T(2,0,1,0)");
  test_for_zero((t4ddg_2(2,0,1,1) - 7*(t1_3(1)*t1_3(1)))/t4ddg_2(2,0,1,1)
		,"T4ddg*T(2,0,1,1)");
  test_for_zero((t4ddg_2(2,0,1,2) - 7*(t1_3(1)*t1_3(2)))/t4ddg_2(2,0,1,2)
		,"T4ddg*T(2,0,1,2)");
  test_for_zero((t4ddg_2(2,0,2,0) - 7*(t1_3(2)*t1_3(0)))/t4ddg_2(2,0,2,0)
		,"T4ddg*T(2,0,2,0)");
  test_for_zero((t4ddg_2(2,0,2,1) - 7*(t1_3(2)*t1_3(1)))/t4ddg_2(2,0,2,1)
		,"T4ddg*T(2,0,2,1)");
  test_for_zero((t4ddg_2(2,0,2,2) - 7*(t1_3(2)*t1_3(2)))/t4ddg_2(2,0,2,2)
		,"T4ddg*T(2,0,2,2)");
  test_for_zero(t4ddg_2(2,1,0,0) - 7*(t1_2(0)*t1_2(0))
		,"T4ddg*T(2,1,0,0)");
  test_for_zero(t4ddg_2(2,1,0,1) - 7*(t1_2(0)*t1_2(1))
		,"T4ddg*T(2,1,0,1)");
  test_for_zero(t4ddg_2(2,1,0,2) - 7*(t1_2(0)*t1_2(2))
		,"T4ddg*T(2,1,0,2)");
  test_for_zero(t4ddg_2(2,1,1,0) - 7*(t1_2(1)*t1_2(0))
		,"T4ddg*T(2,1,1,0)");
  test_for_zero(t4ddg_2(2,1,1,1) - 7*(t1_2(1)*t1_2(1))
		,"T4ddg*T(2,1,1,1)");
  test_for_zero(t4ddg_2(2,1,1,2) - 7*(t1_2(1)*t1_2(2))
		,"T4ddg*T(2,1,1,2)");
  test_for_zero(t4ddg_2(2,1,2,0) - 7*(t1_2(2)*t1_2(0))
		,"T4ddg*T(2,1,2,0)");
  test_for_zero(t4ddg_2(2,1,2,1) - 7*(t1_2(2)*t1_2(1))
		,"T4ddg*T(2,1,2,1)");
  test_for_zero(t4ddg_2(2,1,2,2) - 7*(t1_2(2)*t1_2(2))
		,"T4ddg*T(2,1,2,2)");
  test_for_zero(t4ddg_2(2,2,0,0) - 7*(t2s_3(0,0))
		,"T4ddg*T(2,2,0,0)");
  test_for_zero(t4ddg_2(2,2,0,1) - 7*(t2s_3(0,1))
		,"T4ddg*T(2,2,0,1)");
  test_for_zero(t4ddg_2(2,2,0,2) - 7*(t2s_3(0,2))
		,"T4ddg*T(2,2,0,2)");
  test_for_zero(t4ddg_2(2,2,1,0) - 7*(t2s_3(1,0))
		,"T4ddg*T(2,2,1,0)");
  test_for_zero(t4ddg_2(2,2,1,1) - 7*(t2s_3(1,1))
		,"T4ddg*T(2,2,1,1)");
  test_for_zero(t4ddg_2(2,2,1,2) - 7*(t2s_3(1,2))
		,"T4ddg*T(2,2,1,2)");
  test_for_zero(t4ddg_2(2,2,2,0) - 7*(t2s_3(2,0))
		,"T4ddg*T(2,2,2,0)");
  test_for_zero(t4ddg_2(2,2,2,1) - 7*(t2s_3(2,1))
		,"T4ddg*T(2,2,2,1)");
  test_for_zero(t4ddg_2(2,2,2,2) - 7*(t2s_3(2,2))
		,"T4ddg*T(2,2,2,2)");

  t4ddg_1(N0,i,j,k)=t3dg_1(j,k,i);

  test_for_zero(t4ddg_1(0,0,0,0) - (t3dg_1(0,0,0))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(0,0,0,0)");
  test_for_zero(t4ddg_1(0,0,0,1) - (t3dg_1(0,1,0))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(0,0,0,1)");
  test_for_zero(t4ddg_1(0,0,0,2) - (t3dg_1(0,2,0))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(0,0,0,2)");
  test_for_zero(t4ddg_1(0,0,1,0) - (t3dg_1(1,0,0))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(0,0,1,0)");
  test_for_zero(t4ddg_1(0,0,1,1) - (t3dg_1(1,1,0))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(0,0,1,1)");
  test_for_zero(t4ddg_1(0,0,1,2) - (t3dg_1(1,2,0))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(0,0,1,2)");
  test_for_zero(t4ddg_1(0,0,2,0) - (t3dg_1(2,0,0))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(0,0,2,0)");
  test_for_zero(t4ddg_1(0,0,2,1) - (t3dg_1(2,1,0))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(0,0,2,1)");
  test_for_zero(t4ddg_1(0,0,2,2) - (t3dg_1(2,2,0))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(0,0,2,2)");
  test_for_zero(t4ddg_1(0,1,0,0) - (t3dg_1(0,0,1))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(0,1,0,0)");
  test_for_zero(t4ddg_1(0,1,0,1) - (t3dg_1(0,1,1))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(0,1,0,1)");
  test_for_zero(t4ddg_1(0,1,0,2) - (t3dg_1(0,2,1))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(0,1,0,2)");
  test_for_zero(t4ddg_1(0,1,1,0) - (t3dg_1(1,0,1))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(0,1,1,0)");
  test_for_zero(t4ddg_1(0,1,1,1) - (t3dg_1(1,1,1))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(0,1,1,1)");
  test_for_zero(t4ddg_1(0,1,1,2) - (t3dg_1(1,2,1))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(0,1,1,2)");
  test_for_zero(t4ddg_1(0,1,2,0) - (t3dg_1(2,0,1))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(0,1,2,0)");
  test_for_zero(t4ddg_1(0,1,2,1) - (t3dg_1(2,1,1))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(0,1,2,1)");
  test_for_zero(t4ddg_1(0,1,2,2) - (t3dg_1(2,2,1))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(0,1,2,2)");
  test_for_zero(t4ddg_1(0,2,0,0) - (t3dg_1(0,0,2))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(0,2,0,0)");
  test_for_zero(t4ddg_1(0,2,0,1) - (t3dg_1(0,1,2))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(0,2,0,1)");
  test_for_zero(t4ddg_1(0,2,0,2) - (t3dg_1(0,2,2))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(0,2,0,2)");
  test_for_zero(t4ddg_1(0,2,1,0) - (t3dg_1(1,0,2))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(0,2,1,0)");
  test_for_zero(t4ddg_1(0,2,1,1) - (t3dg_1(1,1,2))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(0,2,1,1)");
  test_for_zero(t4ddg_1(0,2,1,2) - (t3dg_1(1,2,2))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(0,2,1,2)");
  test_for_zero(t4ddg_1(0,2,2,0) - (t3dg_1(2,0,2))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(0,2,2,0)");
  test_for_zero(t4ddg_1(0,2,2,1) - (t3dg_1(2,1,2))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(0,2,2,1)");
  test_for_zero(t4ddg_1(0,2,2,2) - (t3dg_1(2,2,2))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(0,2,2,2)");

  t4ddg_1(N1,i,j,k)=t3dg_2(j,k,i);

  test_for_zero(t4ddg_1(1,0,0,0) - (t3dg_2(0,0,0))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(1,0,0,0)");
  test_for_zero(t4ddg_1(1,0,0,1) - (t3dg_2(0,1,0))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(1,0,0,1)");
  test_for_zero(t4ddg_1(1,0,0,2) - (t3dg_2(0,2,0))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(1,0,0,2)");
  test_for_zero(t4ddg_1(1,0,1,0) - (t3dg_2(1,0,0))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(1,0,1,0)");
  test_for_zero(t4ddg_1(1,0,1,1) - (t3dg_2(1,1,0))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(1,0,1,1)");
  test_for_zero(t4ddg_1(1,0,1,2) - (t3dg_2(1,2,0))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(1,0,1,2)");
  test_for_zero(t4ddg_1(1,0,2,0) - (t3dg_2(2,0,0))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(1,0,2,0)");
  test_for_zero(t4ddg_1(1,0,2,1) - (t3dg_2(2,1,0))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(1,0,2,1)");
  test_for_zero(t4ddg_1(1,0,2,2) - (t3dg_2(2,2,0))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(1,0,2,2)");
  test_for_zero(t4ddg_1(1,1,0,0) - (t3dg_2(0,0,1))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(1,1,0,0)");
  test_for_zero(t4ddg_1(1,1,0,1) - (t3dg_2(0,1,1))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(1,1,0,1)");
  test_for_zero(t4ddg_1(1,1,0,2) - (t3dg_2(0,2,1))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(1,1,0,2)");
  test_for_zero(t4ddg_1(1,1,1,0) - (t3dg_2(1,0,1))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(1,1,1,0)");
  test_for_zero(t4ddg_1(1,1,1,1) - (t3dg_2(1,1,1))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(1,1,1,1)");
  test_for_zero(t4ddg_1(1,1,1,2) - (t3dg_2(1,2,1))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(1,1,1,2)");
  test_for_zero(t4ddg_1(1,1,2,0) - (t3dg_2(2,0,1))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(1,1,2,0)");
  test_for_zero(t4ddg_1(1,1,2,1) - (t3dg_2(2,1,1))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(1,1,2,1)");
  test_for_zero(t4ddg_1(1,1,2,2) - (t3dg_2(2,2,1))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(1,1,2,2)");
  test_for_zero(t4ddg_1(1,2,0,0) - (t3dg_2(0,0,2))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(1,2,0,0)");
  test_for_zero(t4ddg_1(1,2,0,1) - (t3dg_2(0,1,2))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(1,2,0,1)");
  test_for_zero(t4ddg_1(1,2,0,2) - (t3dg_2(0,2,2))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(1,2,0,2)");
  test_for_zero(t4ddg_1(1,2,1,0) - (t3dg_2(1,0,2))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(1,2,1,0)");
  test_for_zero(t4ddg_1(1,2,1,1) - (t3dg_2(1,1,2))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(1,2,1,1)");
  test_for_zero(t4ddg_1(1,2,1,2) - (t3dg_2(1,2,2))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(1,2,1,2)");
  test_for_zero(t4ddg_1(1,2,2,0) - (t3dg_2(2,0,2))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(1,2,2,0)");
  test_for_zero(t4ddg_1(1,2,2,1) - (t3dg_2(2,1,2))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(1,2,2,1)");
  test_for_zero(t4ddg_1(1,2,2,2) - (t3dg_2(2,2,2))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(1,2,2,2)");

  t4ddg_1(N2,i,j,k)=t3dg_3(j,k,i);

  test_for_zero(t4ddg_1(2,0,0,0) - (t3dg_3(0,0,0))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(2,0,0,0)");
  test_for_zero(t4ddg_1(2,0,0,1) - (t3dg_3(0,1,0))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(2,0,0,1)");
  test_for_zero(t4ddg_1(2,0,0,2) - (t3dg_3(0,2,0))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(2,0,0,2)");
  test_for_zero(t4ddg_1(2,0,1,0) - (t3dg_3(1,0,0))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(2,0,1,0)");
  test_for_zero(t4ddg_1(2,0,1,1) - (t3dg_3(1,1,0))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(2,0,1,1)");
  test_for_zero(t4ddg_1(2,0,1,2) - (t3dg_3(1,2,0))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(2,0,1,2)");
  test_for_zero(t4ddg_1(2,0,2,0) - (t3dg_3(2,0,0))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(2,0,2,0)");
  test_for_zero(t4ddg_1(2,0,2,1) - (t3dg_3(2,1,0))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(2,0,2,1)");
  test_for_zero(t4ddg_1(2,0,2,2) - (t3dg_3(2,2,0))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(2,0,2,2)");
  test_for_zero(t4ddg_1(2,1,0,0) - (t3dg_3(0,0,1))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(2,1,0,0)");
  test_for_zero(t4ddg_1(2,1,0,1) - (t3dg_3(0,1,1))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(2,1,0,1)");
  test_for_zero(t4ddg_1(2,1,0,2) - (t3dg_3(0,2,1))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(2,1,0,2)");
  test_for_zero(t4ddg_1(2,1,1,0) - (t3dg_3(1,0,1))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(2,1,1,0)");
  test_for_zero(t4ddg_1(2,1,1,1) - (t3dg_3(1,1,1))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(2,1,1,1)");
  test_for_zero(t4ddg_1(2,1,1,2) - (t3dg_3(1,2,1))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(2,1,1,2)");
  test_for_zero(t4ddg_1(2,1,2,0) - (t3dg_3(2,0,1))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(2,1,2,0)");
  test_for_zero(t4ddg_1(2,1,2,1) - (t3dg_3(2,1,1))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(2,1,2,1)");
  test_for_zero(t4ddg_1(2,1,2,2) - (t3dg_3(2,2,1))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(2,1,2,2)");
  test_for_zero(t4ddg_1(2,2,0,0) - (t3dg_3(0,0,2))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(2,2,0,0)");
  test_for_zero(t4ddg_1(2,2,0,1) - (t3dg_3(0,1,2))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(2,2,0,1)");
  test_for_zero(t4ddg_1(2,2,0,2) - (t3dg_3(0,2,2))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(2,2,0,2)");
  test_for_zero(t4ddg_1(2,2,1,0) - (t3dg_3(1,0,2))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(2,2,1,0)");
  test_for_zero(t4ddg_1(2,2,1,1) - (t3dg_3(1,1,2))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(2,2,1,1)");
  test_for_zero(t4ddg_1(2,2,1,2) - (t3dg_3(1,2,2))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(2,2,1,2)");
  test_for_zero(t4ddg_1(2,2,2,0) - (t3dg_3(2,0,2))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(2,2,2,0)");
  test_for_zero(t4ddg_1(2,2,2,1) - (t3dg_3(2,1,2))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(2,2,2,1)");
  test_for_zero(t4ddg_1(2,2,2,2) - (t3dg_3(2,2,2))
		,"T4ddg(N,i,j,k)=T3dg(j,k,i)(2,2,2,2)");

  t2s_1(i,j)=t4ddg_1(0,0,i,j);
  test_for_zero(t4ddg_1(0,0,0,0) - (t2s_1(0,0))
		,"T4ddg(Num,Num,i,j)(0,0,0,0)");
  test_for_zero(t4ddg_1(0,0,0,1) - (t2s_1(0,1))
		,"T4ddg(Num,Num,i,j)(0,0,0,1)");
  test_for_zero(t4ddg_1(0,0,0,2) - (t2s_1(0,2))
		,"T4ddg(Num,Num,i,j)(0,0,0,2)");
  test_for_zero(t4ddg_1(0,0,1,0) - (t2s_1(1,0))
		,"T4ddg(Num,Num,i,j)(0,0,1,0)");
  test_for_zero(t4ddg_1(0,0,1,1) - (t2s_1(1,1))
		,"T4ddg(Num,Num,i,j)(0,0,1,1)");
  test_for_zero(t4ddg_1(0,0,1,2) - (t2s_1(1,2))
		,"T4ddg(Num,Num,i,j)(0,0,1,2)");
  test_for_zero(t4ddg_1(0,0,2,0) - (t2s_1(2,0))
		,"T4ddg(Num,Num,i,j)(0,0,2,0)");
  test_for_zero(t4ddg_1(0,0,2,1) - (t2s_1(2,1))
		,"T4ddg(Num,Num,i,j)(0,0,2,1)");
  test_for_zero(t4ddg_1(0,0,2,2) - (t2s_1(2,2))
		,"T4ddg(Num,Num,i,j)(0,0,2,2)");

  t2s_1(i,j)=t4ddg_1(0,1,i,j);
  test_for_zero(t4ddg_1(0,1,0,0) - (t2s_1(0,0))
		,"T4ddg(Num,Num,i,j)(0,1,0,0)");
  test_for_zero(t4ddg_1(0,1,0,1) - (t2s_1(0,1))
		,"T4ddg(Num,Num,i,j)(0,1,0,1)");
  test_for_zero(t4ddg_1(0,1,0,2) - (t2s_1(0,2))
		,"T4ddg(Num,Num,i,j)(0,1,0,2)");
  test_for_zero(t4ddg_1(0,1,1,0) - (t2s_1(1,0))
		,"T4ddg(Num,Num,i,j)(0,1,1,0)");
  test_for_zero(t4ddg_1(0,1,1,1) - (t2s_1(1,1))
		,"T4ddg(Num,Num,i,j)(0,1,1,1)");
  test_for_zero(t4ddg_1(0,1,1,2) - (t2s_1(1,2))
		,"T4ddg(Num,Num,i,j)(0,1,1,2)");
  test_for_zero(t4ddg_1(0,1,2,0) - (t2s_1(2,0))
		,"T4ddg(Num,Num,i,j)(0,1,2,0)");
  test_for_zero(t4ddg_1(0,1,2,1) - (t2s_1(2,1))
		,"T4ddg(Num,Num,i,j)(0,1,2,1)");
  test_for_zero(t4ddg_1(0,1,2,2) - (t2s_1(2,2))
		,"T4ddg(Num,Num,i,j)(0,1,2,2)");

  t2s_1(i,j)=t4ddg_1(0,2,i,j);
  test_for_zero(t4ddg_1(0,2,0,0) - (t2s_1(0,0))
		,"T4ddg(Num,Num,i,j)(0,2,0,0)");
  test_for_zero(t4ddg_1(0,2,0,1) - (t2s_1(0,1))
		,"T4ddg(Num,Num,i,j)(0,2,0,1)");
  test_for_zero(t4ddg_1(0,2,0,2) - (t2s_1(0,2))
		,"T4ddg(Num,Num,i,j)(0,2,0,2)");
  test_for_zero(t4ddg_1(0,2,1,0) - (t2s_1(1,0))
		,"T4ddg(Num,Num,i,j)(0,2,1,0)");
  test_for_zero(t4ddg_1(0,2,1,1) - (t2s_1(1,1))
		,"T4ddg(Num,Num,i,j)(0,2,1,1)");
  test_for_zero(t4ddg_1(0,2,1,2) - (t2s_1(1,2))
		,"T4ddg(Num,Num,i,j)(0,2,1,2)");
  test_for_zero(t4ddg_1(0,2,2,0) - (t2s_1(2,0))
		,"T4ddg(Num,Num,i,j)(0,2,2,0)");
  test_for_zero(t4ddg_1(0,2,2,1) - (t2s_1(2,1))
		,"T4ddg(Num,Num,i,j)(0,2,2,1)");
  test_for_zero(t4ddg_1(0,2,2,2) - (t2s_1(2,2))
		,"T4ddg(Num,Num,i,j)(0,2,2,2)");

  t2s_1(i,j)=t4ddg_1(1,0,i,j);
  test_for_zero(t4ddg_1(1,0,0,0) - (t2s_1(0,0))
		,"T4ddg(Num,Num,i,j)(1,0,0,0)");
  test_for_zero(t4ddg_1(1,0,0,1) - (t2s_1(0,1))
		,"T4ddg(Num,Num,i,j)(1,0,0,1)");
  test_for_zero(t4ddg_1(1,0,0,2) - (t2s_1(0,2))
		,"T4ddg(Num,Num,i,j)(1,0,0,2)");
  test_for_zero(t4ddg_1(1,0,1,0) - (t2s_1(1,0))
		,"T4ddg(Num,Num,i,j)(1,0,1,0)");
  test_for_zero(t4ddg_1(1,0,1,1) - (t2s_1(1,1))
		,"T4ddg(Num,Num,i,j)(1,0,1,1)");
  test_for_zero(t4ddg_1(1,0,1,2) - (t2s_1(1,2))
		,"T4ddg(Num,Num,i,j)(1,0,1,2)");
  test_for_zero(t4ddg_1(1,0,2,0) - (t2s_1(2,0))
		,"T4ddg(Num,Num,i,j)(1,0,2,0)");
  test_for_zero(t4ddg_1(1,0,2,1) - (t2s_1(2,1))
		,"T4ddg(Num,Num,i,j)(1,0,2,1)");
  test_for_zero(t4ddg_1(1,0,2,2) - (t2s_1(2,2))
		,"T4ddg(Num,Num,i,j)(1,0,2,2)");

  t2s_1(i,j)=t4ddg_1(1,1,i,j);
  test_for_zero(t4ddg_1(1,1,0,0) - (t2s_1(0,0))
		,"T4ddg(Num,Num,i,j)(1,1,0,0)");
  test_for_zero(t4ddg_1(1,1,0,1) - (t2s_1(0,1))
		,"T4ddg(Num,Num,i,j)(1,1,0,1)");
  test_for_zero(t4ddg_1(1,1,0,2) - (t2s_1(0,2))
		,"T4ddg(Num,Num,i,j)(1,1,0,2)");
  test_for_zero(t4ddg_1(1,1,1,0) - (t2s_1(1,0))
		,"T4ddg(Num,Num,i,j)(1,1,1,0)");
  test_for_zero(t4ddg_1(1,1,1,1) - (t2s_1(1,1))
		,"T4ddg(Num,Num,i,j)(1,1,1,1)");
  test_for_zero(t4ddg_1(1,1,1,2) - (t2s_1(1,2))
		,"T4ddg(Num,Num,i,j)(1,1,1,2)");
  test_for_zero(t4ddg_1(1,1,2,0) - (t2s_1(2,0))
		,"T4ddg(Num,Num,i,j)(1,1,2,0)");
  test_for_zero(t4ddg_1(1,1,2,1) - (t2s_1(2,1))
		,"T4ddg(Num,Num,i,j)(1,1,2,1)");
  test_for_zero(t4ddg_1(1,1,2,2) - (t2s_1(2,2))
		,"T4ddg(Num,Num,i,j)(1,1,2,2)");

  t2s_1(i,j)=t4ddg_1(1,2,i,j);
  test_for_zero(t4ddg_1(1,2,0,0) - (t2s_1(0,0))
		,"T4ddg(Num,Num,i,j)(1,2,0,0)");
  test_for_zero(t4ddg_1(1,2,0,1) - (t2s_1(0,1))
		,"T4ddg(Num,Num,i,j)(1,2,0,1)");
  test_for_zero(t4ddg_1(1,2,0,2) - (t2s_1(0,2))
		,"T4ddg(Num,Num,i,j)(1,2,0,2)");
  test_for_zero(t4ddg_1(1,2,1,0) - (t2s_1(1,0))
		,"T4ddg(Num,Num,i,j)(1,2,1,0)");
  test_for_zero(t4ddg_1(1,2,1,1) - (t2s_1(1,1))
		,"T4ddg(Num,Num,i,j)(1,2,1,1)");
  test_for_zero(t4ddg_1(1,2,1,2) - (t2s_1(1,2))
		,"T4ddg(Num,Num,i,j)(1,2,1,2)");
  test_for_zero(t4ddg_1(1,2,2,0) - (t2s_1(2,0))
		,"T4ddg(Num,Num,i,j)(1,2,2,0)");
  test_for_zero(t4ddg_1(1,2,2,1) - (t2s_1(2,1))
		,"T4ddg(Num,Num,i,j)(1,2,2,1)");
  test_for_zero(t4ddg_1(1,2,2,2) - (t2s_1(2,2))
		,"T4ddg(Num,Num,i,j)(1,2,2,2)");

  t2s_1(i,j)=t4ddg_1(2,0,i,j);
  test_for_zero(t4ddg_1(2,0,0,0) - (t2s_1(0,0))
		,"T4ddg(Num,Num,i,j)(2,0,0,0)");
  test_for_zero(t4ddg_1(2,0,0,1) - (t2s_1(0,1))
		,"T4ddg(Num,Num,i,j)(2,0,0,1)");
  test_for_zero(t4ddg_1(2,0,0,2) - (t2s_1(0,2))
		,"T4ddg(Num,Num,i,j)(2,0,0,2)");
  test_for_zero(t4ddg_1(2,0,1,0) - (t2s_1(1,0))
		,"T4ddg(Num,Num,i,j)(2,0,1,0)");
  test_for_zero(t4ddg_1(2,0,1,1) - (t2s_1(1,1))
		,"T4ddg(Num,Num,i,j)(2,0,1,1)");
  test_for_zero(t4ddg_1(2,0,1,2) - (t2s_1(1,2))
		,"T4ddg(Num,Num,i,j)(2,0,1,2)");
  test_for_zero(t4ddg_1(2,0,2,0) - (t2s_1(2,0))
		,"T4ddg(Num,Num,i,j)(2,0,2,0)");
  test_for_zero(t4ddg_1(2,0,2,1) - (t2s_1(2,1))
		,"T4ddg(Num,Num,i,j)(2,0,2,1)");
  test_for_zero(t4ddg_1(2,0,2,2) - (t2s_1(2,2))
		,"T4ddg(Num,Num,i,j)(2,0,2,2)");

  t2s_1(i,j)=t4ddg_1(2,1,i,j);
  test_for_zero(t4ddg_1(2,1,0,0) - (t2s_1(0,0))
		,"T4ddg(Num,Num,i,j)(2,1,0,0)");
  test_for_zero(t4ddg_1(2,1,0,1) - (t2s_1(0,1))
		,"T4ddg(Num,Num,i,j)(2,1,0,1)");
  test_for_zero(t4ddg_1(2,1,0,2) - (t2s_1(0,2))
		,"T4ddg(Num,Num,i,j)(2,1,0,2)");
  test_for_zero(t4ddg_1(2,1,1,0) - (t2s_1(1,0))
		,"T4ddg(Num,Num,i,j)(2,1,1,0)");
  test_for_zero(t4ddg_1(2,1,1,1) - (t2s_1(1,1))
		,"T4ddg(Num,Num,i,j)(2,1,1,1)");
  test_for_zero(t4ddg_1(2,1,1,2) - (t2s_1(1,2))
		,"T4ddg(Num,Num,i,j)(2,1,1,2)");
  test_for_zero(t4ddg_1(2,1,2,0) - (t2s_1(2,0))
		,"T4ddg(Num,Num,i,j)(2,1,2,0)");
  test_for_zero(t4ddg_1(2,1,2,1) - (t2s_1(2,1))
		,"T4ddg(Num,Num,i,j)(2,1,2,1)");
  test_for_zero(t4ddg_1(2,1,2,2) - (t2s_1(2,2))
		,"T4ddg(Num,Num,i,j)(2,1,2,2)");

  t2s_1(i,j)=t4ddg_1(2,2,i,j);
  test_for_zero(t4ddg_1(2,2,0,0) - (t2s_1(0,0))
		,"T4ddg(Num,Num,i,j)(2,2,0,0)");
  test_for_zero(t4ddg_1(2,2,0,1) - (t2s_1(0,1))
		,"T4ddg(Num,Num,i,j)(2,2,0,1)");
  test_for_zero(t4ddg_1(2,2,0,2) - (t2s_1(0,2))
		,"T4ddg(Num,Num,i,j)(2,2,0,2)");
  test_for_zero(t4ddg_1(2,2,1,0) - (t2s_1(1,0))
		,"T4ddg(Num,Num,i,j)(2,2,1,0)");
  test_for_zero(t4ddg_1(2,2,1,1) - (t2s_1(1,1))
		,"T4ddg(Num,Num,i,j)(2,2,1,1)");
  test_for_zero(t4ddg_1(2,2,1,2) - (t2s_1(1,2))
		,"T4ddg(Num,Num,i,j)(2,2,1,2)");
  test_for_zero(t4ddg_1(2,2,2,0) - (t2s_1(2,0))
		,"T4ddg(Num,Num,i,j)(2,2,2,0)");
  test_for_zero(t4ddg_1(2,2,2,1) - (t2s_1(2,1))
		,"T4ddg(Num,Num,i,j)(2,2,2,1)");
  test_for_zero(t4ddg_1(2,2,2,2) - (t2s_1(2,2))
		,"T4ddg(Num,Num,i,j)(2,2,2,2)");





  t2s_1(i,j)=t4ddg_1(i,j,0,0);
  test_for_zero(t4ddg_1(0,0,0,0) - (t2s_1(0,0))
		,"T4ddg(i,j,Num,Num)(0,0,0,0)");
  test_for_zero(t4ddg_1(0,1,0,0) - (t2s_1(0,1))
		,"T4ddg(i,j,Num,Num)(0,0,0,1)");
  test_for_zero(t4ddg_1(0,2,0,0) - (t2s_1(0,2))
		,"T4ddg(i,j,Num,Num)(0,0,0,2)");
  test_for_zero(t4ddg_1(1,0,0,0) - (t2s_1(1,0))
		,"T4ddg(i,j,Num,Num)(0,0,1,0)");
  test_for_zero(t4ddg_1(1,1,0,0) - (t2s_1(1,1))
		,"T4ddg(i,j,Num,Num)(0,0,1,1)");
  test_for_zero(t4ddg_1(1,2,0,0) - (t2s_1(1,2))
		,"T4ddg(i,j,Num,Num)(0,0,1,2)");
  test_for_zero(t4ddg_1(2,0,0,0) - (t2s_1(2,0))
		,"T4ddg(i,j,Num,Num)(0,0,2,0)");
  test_for_zero(t4ddg_1(2,1,0,0) - (t2s_1(2,1))
		,"T4ddg(i,j,Num,Num)(0,0,2,1)");
  test_for_zero(t4ddg_1(2,2,0,0) - (t2s_1(2,2))
		,"T4ddg(i,j,Num,Num)(0,0,2,2)");

  t2s_1(i,j)=t4ddg_1(i,j,0,1);
  test_for_zero(t4ddg_1(0,0,0,1) - (t2s_1(0,0))
		,"T4ddg(i,j,Num,Num)(0,1,0,0)");
  test_for_zero(t4ddg_1(0,1,0,1) - (t2s_1(0,1))
		,"T4ddg(i,j,Num,Num)(0,1,0,1)");
  test_for_zero(t4ddg_1(0,2,0,1) - (t2s_1(0,2))
		,"T4ddg(i,j,Num,Num)(0,1,0,2)");
  test_for_zero(t4ddg_1(1,0,0,1) - (t2s_1(1,0))
		,"T4ddg(i,j,Num,Num)(0,1,1,0)");
  test_for_zero(t4ddg_1(1,1,0,1) - (t2s_1(1,1))
		,"T4ddg(i,j,Num,Num)(0,1,1,1)");
  test_for_zero(t4ddg_1(1,2,0,1) - (t2s_1(1,2))
		,"T4ddg(i,j,Num,Num)(0,1,1,2)");
  test_for_zero(t4ddg_1(2,0,0,1) - (t2s_1(2,0))
		,"T4ddg(i,j,Num,Num)(0,1,2,0)");
  test_for_zero(t4ddg_1(2,1,0,1) - (t2s_1(2,1))
		,"T4ddg(i,j,Num,Num)(0,1,2,1)");
  test_for_zero(t4ddg_1(2,2,0,1) - (t2s_1(2,2))
		,"T4ddg(i,j,Num,Num)(0,1,2,2)");

  t2s_1(i,j)=t4ddg_1(i,j,0,2);
  test_for_zero(t4ddg_1(0,0,0,2) - (t2s_1(0,0))
		,"T4ddg(i,j,Num,Num)(0,2,0,0)");
  test_for_zero(t4ddg_1(0,1,0,2) - (t2s_1(0,1))
		,"T4ddg(i,j,Num,Num)(0,2,0,1)");
  test_for_zero(t4ddg_1(0,2,0,2) - (t2s_1(0,2))
		,"T4ddg(i,j,Num,Num)(0,2,0,2)");
  test_for_zero(t4ddg_1(1,0,0,2) - (t2s_1(1,0))
		,"T4ddg(i,j,Num,Num)(0,2,1,0)");
  test_for_zero(t4ddg_1(1,1,0,2) - (t2s_1(1,1))
		,"T4ddg(i,j,Num,Num)(0,2,1,1)");
  test_for_zero(t4ddg_1(1,2,0,2) - (t2s_1(1,2))
		,"T4ddg(i,j,Num,Num)(0,2,1,2)");
  test_for_zero(t4ddg_1(2,0,0,2) - (t2s_1(2,0))
		,"T4ddg(i,j,Num,Num)(0,2,2,0)");
  test_for_zero(t4ddg_1(2,1,0,2) - (t2s_1(2,1))
		,"T4ddg(i,j,Num,Num)(0,2,2,1)");
  test_for_zero(t4ddg_1(2,2,0,2) - (t2s_1(2,2))
		,"T4ddg(i,j,Num,Num)(0,2,2,2)");

  t2s_1(i,j)=t4ddg_1(i,j,1,0);
  test_for_zero(t4ddg_1(0,0,1,0) - (t2s_1(0,0))
		,"T4ddg(i,j,Num,Num)(1,0,0,0)");
  test_for_zero(t4ddg_1(0,1,1,0) - (t2s_1(0,1))
		,"T4ddg(i,j,Num,Num)(1,0,0,1)");
  test_for_zero(t4ddg_1(0,2,1,0) - (t2s_1(0,2))
		,"T4ddg(i,j,Num,Num)(1,0,0,2)");
  test_for_zero(t4ddg_1(1,0,1,0) - (t2s_1(1,0))
		,"T4ddg(i,j,Num,Num)(1,0,1,0)");
  test_for_zero(t4ddg_1(1,1,1,0) - (t2s_1(1,1))
		,"T4ddg(i,j,Num,Num)(1,0,1,1)");
  test_for_zero(t4ddg_1(1,2,1,0) - (t2s_1(1,2))
		,"T4ddg(i,j,Num,Num)(1,0,1,2)");
  test_for_zero(t4ddg_1(2,0,1,0) - (t2s_1(2,0))
		,"T4ddg(i,j,Num,Num)(1,0,2,0)");
  test_for_zero(t4ddg_1(2,1,1,0) - (t2s_1(2,1))
		,"T4ddg(i,j,Num,Num)(1,0,2,1)");
  test_for_zero(t4ddg_1(2,2,1,0) - (t2s_1(2,2))
		,"T4ddg(i,j,Num,Num)(1,0,2,2)");

  t2s_1(i,j)=t4ddg_1(i,j,1,1);
  test_for_zero(t4ddg_1(0,0,1,1) - (t2s_1(0,0))
		,"T4ddg(i,j,Num,Num)(1,1,0,0)");
  test_for_zero(t4ddg_1(0,1,1,1) - (t2s_1(0,1))
		,"T4ddg(i,j,Num,Num)(1,1,0,1)");
  test_for_zero(t4ddg_1(0,2,1,1) - (t2s_1(0,2))
		,"T4ddg(i,j,Num,Num)(1,1,0,2)");
  test_for_zero(t4ddg_1(1,0,1,1) - (t2s_1(1,0))
		,"T4ddg(i,j,Num,Num)(1,1,1,0)");
  test_for_zero(t4ddg_1(1,1,1,1) - (t2s_1(1,1))
		,"T4ddg(i,j,Num,Num)(1,1,1,1)");
  test_for_zero(t4ddg_1(1,2,1,1) - (t2s_1(1,2))
		,"T4ddg(i,j,Num,Num)(1,1,1,2)");
  test_for_zero(t4ddg_1(2,0,1,1) - (t2s_1(2,0))
		,"T4ddg(i,j,Num,Num)(1,1,2,0)");
  test_for_zero(t4ddg_1(2,1,1,1) - (t2s_1(2,1))
		,"T4ddg(i,j,Num,Num)(1,1,2,1)");
  test_for_zero(t4ddg_1(2,2,1,1) - (t2s_1(2,2))
		,"T4ddg(i,j,Num,Num)(1,1,2,2)");

  t2s_1(i,j)=t4ddg_1(i,j,1,2);
  test_for_zero(t4ddg_1(0,0,1,2) - (t2s_1(0,0))
		,"T4ddg(i,j,Num,Num)(1,2,0,0)");
  test_for_zero(t4ddg_1(0,1,1,2) - (t2s_1(0,1))
		,"T4ddg(i,j,Num,Num)(1,2,0,1)");
  test_for_zero(t4ddg_1(0,2,1,2) - (t2s_1(0,2))
		,"T4ddg(i,j,Num,Num)(1,2,0,2)");
  test_for_zero(t4ddg_1(1,0,1,2) - (t2s_1(1,0))
		,"T4ddg(i,j,Num,Num)(1,2,1,0)");
  test_for_zero(t4ddg_1(1,1,1,2) - (t2s_1(1,1))
		,"T4ddg(i,j,Num,Num)(1,2,1,1)");
  test_for_zero(t4ddg_1(1,2,1,2) - (t2s_1(1,2))
		,"T4ddg(i,j,Num,Num)(1,2,1,2)");
  test_for_zero(t4ddg_1(2,0,1,2) - (t2s_1(2,0))
		,"T4ddg(i,j,Num,Num)(1,2,2,0)");
  test_for_zero(t4ddg_1(2,1,1,2) - (t2s_1(2,1))
		,"T4ddg(i,j,Num,Num)(1,2,2,1)");
  test_for_zero(t4ddg_1(2,2,1,2) - (t2s_1(2,2))
		,"T4ddg(i,j,Num,Num)(1,2,2,2)");

  t2s_1(i,j)=t4ddg_1(i,j,2,0);
  test_for_zero(t4ddg_1(0,0,2,0) - (t2s_1(0,0))
		,"T4ddg(i,j,Num,Num)(2,0,0,0)");
  test_for_zero(t4ddg_1(0,1,2,0) - (t2s_1(0,1))
		,"T4ddg(i,j,Num,Num)(2,0,0,1)");
  test_for_zero(t4ddg_1(0,2,2,0) - (t2s_1(0,2))
		,"T4ddg(i,j,Num,Num)(2,0,0,2)");
  test_for_zero(t4ddg_1(1,0,2,0) - (t2s_1(1,0))
		,"T4ddg(i,j,Num,Num)(2,0,1,0)");
  test_for_zero(t4ddg_1(1,1,2,0) - (t2s_1(1,1))
		,"T4ddg(i,j,Num,Num)(2,0,1,1)");
  test_for_zero(t4ddg_1(1,2,2,0) - (t2s_1(1,2))
		,"T4ddg(i,j,Num,Num)(2,0,1,2)");
  test_for_zero(t4ddg_1(2,0,2,0) - (t2s_1(2,0))
		,"T4ddg(i,j,Num,Num)(2,0,2,0)");
  test_for_zero(t4ddg_1(2,1,2,0) - (t2s_1(2,1))
		,"T4ddg(i,j,Num,Num)(2,0,2,1)");
  test_for_zero(t4ddg_1(2,2,2,0) - (t2s_1(2,2))
		,"T4ddg(i,j,Num,Num)(2,0,2,2)");

  t2s_1(i,j)=t4ddg_1(i,j,2,1);
  test_for_zero(t4ddg_1(0,0,2,1) - (t2s_1(0,0))
		,"T4ddg(i,j,Num,Num)(2,1,0,0)");
  test_for_zero(t4ddg_1(0,1,2,1) - (t2s_1(0,1))
		,"T4ddg(i,j,Num,Num)(2,1,0,1)");
  test_for_zero(t4ddg_1(0,2,2,1) - (t2s_1(0,2))
		,"T4ddg(i,j,Num,Num)(2,1,0,2)");
  test_for_zero(t4ddg_1(1,0,2,1) - (t2s_1(1,0))
		,"T4ddg(i,j,Num,Num)(2,1,1,0)");
  test_for_zero(t4ddg_1(1,1,2,1) - (t2s_1(1,1))
		,"T4ddg(i,j,Num,Num)(2,1,1,1)");
  test_for_zero(t4ddg_1(1,2,2,1) - (t2s_1(1,2))
		,"T4ddg(i,j,Num,Num)(2,1,1,2)");
  test_for_zero(t4ddg_1(2,0,2,1) - (t2s_1(2,0))
		,"T4ddg(i,j,Num,Num)(2,1,2,0)");
  test_for_zero(t4ddg_1(2,1,2,1) - (t2s_1(2,1))
		,"T4ddg(i,j,Num,Num)(2,1,2,1)");
  test_for_zero(t4ddg_1(2,2,2,1) - (t2s_1(2,2))
		,"T4ddg(i,j,Num,Num)(2,1,2,2)");

  t2s_1(i,j)=t4ddg_1(i,j,2,2);
  test_for_zero(t4ddg_1(0,0,2,2) - (t2s_1(0,0))
		,"T4ddg(i,j,Num,Num)(2,2,0,0)");
  test_for_zero(t4ddg_1(0,1,2,2) - (t2s_1(0,1))
		,"T4ddg(i,j,Num,Num)(2,2,0,1)");
  test_for_zero(t4ddg_1(0,2,2,2) - (t2s_1(0,2))
		,"T4ddg(i,j,Num,Num)(2,2,0,2)");
  test_for_zero(t4ddg_1(1,0,2,2) - (t2s_1(1,0))
		,"T4ddg(i,j,Num,Num)(2,2,1,0)");
  test_for_zero(t4ddg_1(1,1,2,2) - (t2s_1(1,1))
		,"T4ddg(i,j,Num,Num)(2,2,1,1)");
  test_for_zero(t4ddg_1(1,2,2,2) - (t2s_1(1,2))
		,"T4ddg(i,j,Num,Num)(2,2,1,2)");
  test_for_zero(t4ddg_1(2,0,2,2) - (t2s_1(2,0))
		,"T4ddg(i,j,Num,Num)(2,2,2,0)");
  test_for_zero(t4ddg_1(2,1,2,2) - (t2s_1(2,1))
		,"T4ddg(i,j,Num,Num)(2,2,2,1)");
  test_for_zero(t4ddg_1(2,2,2,2) - (t2s_1(2,2))
		,"T4ddg(i,j,Num,Num)(2,2,2,2)");



  t2_1(i,j)=t4ddg_1(0,i,0,j);
  test_for_zero(t4ddg_1(0,0,0,0) - (t2_1(0,0))
		,"T4ddg(Num,i,Num,j)(0,0,0,0)");
  test_for_zero(t4ddg_1(0,0,0,1) - (t2_1(0,1))
		,"T4ddg(Num,i,Num,j)(0,0,0,1)");
  test_for_zero(t4ddg_1(0,0,0,2) - (t2_1(0,2))
		,"T4ddg(Num,i,Num,j)(0,0,0,2)");
  test_for_zero(t4ddg_1(0,1,0,0) - (t2_1(1,0))
		,"T4ddg(Num,i,Num,j)(0,0,1,0)");
  test_for_zero(t4ddg_1(0,1,0,1) - (t2_1(1,1))
		,"T4ddg(Num,i,Num,j)(0,0,1,1)");
  test_for_zero(t4ddg_1(0,1,0,2) - (t2_1(1,2))
		,"T4ddg(Num,i,Num,j)(0,0,1,2)");
  test_for_zero(t4ddg_1(0,2,0,0) - (t2_1(2,0))
		,"T4ddg(Num,i,Num,j)(0,0,2,0)");
  test_for_zero(t4ddg_1(0,2,0,1) - (t2_1(2,1))
		,"T4ddg(Num,i,Num,j)(0,0,2,1)");
  test_for_zero(t4ddg_1(0,2,0,2) - (t2_1(2,2))
		,"T4ddg(Num,i,Num,j)(0,0,2,2)");

  t2_1(i,j)=t4ddg_1(0,i,1,j);
  test_for_zero(t4ddg_1(0,0,1,0) - (t2_1(0,0))
		,"T4ddg(Num,i,Num,j)(0,1,0,0)");
  test_for_zero(t4ddg_1(0,0,1,1) - (t2_1(0,1))
		,"T4ddg(Num,i,Num,j)(0,1,0,1)");
  test_for_zero(t4ddg_1(0,0,1,2) - (t2_1(0,2))
		,"T4ddg(Num,i,Num,j)(0,1,0,2)");
  test_for_zero(t4ddg_1(0,1,1,0) - (t2_1(1,0))
		,"T4ddg(Num,i,Num,j)(0,1,1,0)");
  test_for_zero(t4ddg_1(0,1,1,1) - (t2_1(1,1))
		,"T4ddg(Num,i,Num,j)(0,1,1,1)");
  test_for_zero(t4ddg_1(0,1,1,2) - (t2_1(1,2))
		,"T4ddg(Num,i,Num,j)(0,1,1,2)");
  test_for_zero(t4ddg_1(0,2,1,0) - (t2_1(2,0))
		,"T4ddg(Num,i,Num,j)(0,1,2,0)");
  test_for_zero(t4ddg_1(0,2,1,1) - (t2_1(2,1))
		,"T4ddg(Num,i,Num,j)(0,1,2,1)");
  test_for_zero(t4ddg_1(0,2,1,2) - (t2_1(2,2))
		,"T4ddg(Num,i,Num,j)(0,1,2,2)");

  t2_1(i,j)=t4ddg_1(0,i,2,j);
  test_for_zero(t4ddg_1(0,0,2,0) - (t2_1(0,0))
		,"T4ddg(Num,i,Num,j)(0,2,0,0)");
  test_for_zero(t4ddg_1(0,0,2,1) - (t2_1(0,1))
		,"T4ddg(Num,i,Num,j)(0,2,0,1)");
  test_for_zero(t4ddg_1(0,0,2,2) - (t2_1(0,2))
		,"T4ddg(Num,i,Num,j)(0,2,0,2)");
  test_for_zero(t4ddg_1(0,1,2,0) - (t2_1(1,0))
		,"T4ddg(Num,i,Num,j)(0,2,1,0)");
  test_for_zero(t4ddg_1(0,1,2,1) - (t2_1(1,1))
		,"T4ddg(Num,i,Num,j)(0,2,1,1)");
  test_for_zero(t4ddg_1(0,1,2,2) - (t2_1(1,2))
		,"T4ddg(Num,i,Num,j)(0,2,1,2)");
  test_for_zero(t4ddg_1(0,2,2,0) - (t2_1(2,0))
		,"T4ddg(Num,i,Num,j)(0,2,2,0)");
  test_for_zero(t4ddg_1(0,2,2,1) - (t2_1(2,1))
		,"T4ddg(Num,i,Num,j)(0,2,2,1)");
  test_for_zero(t4ddg_1(0,2,2,2) - (t2_1(2,2))
		,"T4ddg(Num,i,Num,j)(0,2,2,2)");

  t2_1(i,j)=t4ddg_1(1,i,0,j);
  test_for_zero(t4ddg_1(1,0,0,0) - (t2_1(0,0))
		,"T4ddg(Num,i,Num,j)(1,0,0,0)");
  test_for_zero(t4ddg_1(1,0,0,1) - (t2_1(0,1))
		,"T4ddg(Num,i,Num,j)(1,0,0,1)");
  test_for_zero(t4ddg_1(1,0,0,2) - (t2_1(0,2))
		,"T4ddg(Num,i,Num,j)(1,0,0,2)");
  test_for_zero(t4ddg_1(1,1,0,0) - (t2_1(1,0))
		,"T4ddg(Num,i,Num,j)(1,0,1,0)");
  test_for_zero(t4ddg_1(1,1,0,1) - (t2_1(1,1))
		,"T4ddg(Num,i,Num,j)(1,0,1,1)");
  test_for_zero(t4ddg_1(1,1,0,2) - (t2_1(1,2))
		,"T4ddg(Num,i,Num,j)(1,0,1,2)");
  test_for_zero(t4ddg_1(1,2,0,0) - (t2_1(2,0))
		,"T4ddg(Num,i,Num,j)(1,0,2,0)");
  test_for_zero(t4ddg_1(1,2,0,1) - (t2_1(2,1))
		,"T4ddg(Num,i,Num,j)(1,0,2,1)");
  test_for_zero(t4ddg_1(1,2,0,2) - (t2_1(2,2))
		,"T4ddg(Num,i,Num,j)(1,0,2,2)");

  t2_1(i,j)=t4ddg_1(1,i,1,j);
  test_for_zero(t4ddg_1(1,0,1,0) - (t2_1(0,0))
		,"T4ddg(Num,i,Num,j)(1,1,0,0)");
  test_for_zero(t4ddg_1(1,0,1,1) - (t2_1(0,1))
		,"T4ddg(Num,i,Num,j)(1,1,0,1)");
  test_for_zero(t4ddg_1(1,0,1,2) - (t2_1(0,2))
		,"T4ddg(Num,i,Num,j)(1,1,0,2)");
  test_for_zero(t4ddg_1(1,1,1,0) - (t2_1(1,0))
		,"T4ddg(Num,i,Num,j)(1,1,1,0)");
  test_for_zero(t4ddg_1(1,1,1,1) - (t2_1(1,1))
		,"T4ddg(Num,i,Num,j)(1,1,1,1)");
  test_for_zero(t4ddg_1(1,1,1,2) - (t2_1(1,2))
		,"T4ddg(Num,i,Num,j)(1,1,1,2)");
  test_for_zero(t4ddg_1(1,2,1,0) - (t2_1(2,0))
		,"T4ddg(Num,i,Num,j)(1,1,2,0)");
  test_for_zero(t4ddg_1(1,2,1,1) - (t2_1(2,1))
		,"T4ddg(Num,i,Num,j)(1,1,2,1)");
  test_for_zero(t4ddg_1(1,2,1,2) - (t2_1(2,2))
		,"T4ddg(Num,i,Num,j)(1,1,2,2)");

  t2_1(i,j)=t4ddg_1(1,i,2,j);
  test_for_zero(t4ddg_1(1,0,2,0) - (t2_1(0,0))
		,"T4ddg(Num,i,Num,j)(1,2,0,0)");
  test_for_zero(t4ddg_1(1,0,2,1) - (t2_1(0,1))
		,"T4ddg(Num,i,Num,j)(1,2,0,1)");
  test_for_zero(t4ddg_1(1,0,2,2) - (t2_1(0,2))
		,"T4ddg(Num,i,Num,j)(1,2,0,2)");
  test_for_zero(t4ddg_1(1,1,2,0) - (t2_1(1,0))
		,"T4ddg(Num,i,Num,j)(1,2,1,0)");
  test_for_zero(t4ddg_1(1,1,2,1) - (t2_1(1,1))
		,"T4ddg(Num,i,Num,j)(1,2,1,1)");
  test_for_zero(t4ddg_1(1,1,2,2) - (t2_1(1,2))
		,"T4ddg(Num,i,Num,j)(1,2,1,2)");
  test_for_zero(t4ddg_1(1,2,2,0) - (t2_1(2,0))
		,"T4ddg(Num,i,Num,j)(1,2,2,0)");
  test_for_zero(t4ddg_1(1,2,2,1) - (t2_1(2,1))
		,"T4ddg(Num,i,Num,j)(1,2,2,1)");
  test_for_zero(t4ddg_1(1,2,2,2) - (t2_1(2,2))
		,"T4ddg(Num,i,Num,j)(1,2,2,2)");

  t2_1(i,j)=t4ddg_1(2,i,0,j);
  test_for_zero(t4ddg_1(2,0,0,0) - (t2_1(0,0))
		,"T4ddg(Num,i,Num,j)(2,0,0,0)");
  test_for_zero(t4ddg_1(2,0,0,1) - (t2_1(0,1))
		,"T4ddg(Num,i,Num,j)(2,0,0,1)");
  test_for_zero(t4ddg_1(2,0,0,2) - (t2_1(0,2))
		,"T4ddg(Num,i,Num,j)(2,0,0,2)");
  test_for_zero(t4ddg_1(2,1,0,0) - (t2_1(1,0))
		,"T4ddg(Num,i,Num,j)(2,0,1,0)");
  test_for_zero(t4ddg_1(2,1,0,1) - (t2_1(1,1))
		,"T4ddg(Num,i,Num,j)(2,0,1,1)");
  test_for_zero(t4ddg_1(2,1,0,2) - (t2_1(1,2))
		,"T4ddg(Num,i,Num,j)(2,0,1,2)");
  test_for_zero(t4ddg_1(2,2,0,0) - (t2_1(2,0))
		,"T4ddg(Num,i,Num,j)(2,0,2,0)");
  test_for_zero(t4ddg_1(2,2,0,1) - (t2_1(2,1))
		,"T4ddg(Num,i,Num,j)(2,0,2,1)");
  test_for_zero(t4ddg_1(2,2,0,2) - (t2_1(2,2))
		,"T4ddg(Num,i,Num,j)(2,0,2,2)");

  t2_1(i,j)=t4ddg_1(2,i,1,j);
  test_for_zero(t4ddg_1(2,0,1,0) - (t2_1(0,0))
		,"T4ddg(Num,i,Num,j)(2,1,0,0)");
  test_for_zero(t4ddg_1(2,0,1,1) - (t2_1(0,1))
		,"T4ddg(Num,i,Num,j)(2,1,0,1)");
  test_for_zero(t4ddg_1(2,0,1,2) - (t2_1(0,2))
		,"T4ddg(Num,i,Num,j)(2,1,0,2)");
  test_for_zero(t4ddg_1(2,1,1,0) - (t2_1(1,0))
		,"T4ddg(Num,i,Num,j)(2,1,1,0)");
  test_for_zero(t4ddg_1(2,1,1,1) - (t2_1(1,1))
		,"T4ddg(Num,i,Num,j)(2,1,1,1)");
  test_for_zero(t4ddg_1(2,1,1,2) - (t2_1(1,2))
		,"T4ddg(Num,i,Num,j)(2,1,1,2)");
  test_for_zero(t4ddg_1(2,2,1,0) - (t2_1(2,0))
		,"T4ddg(Num,i,Num,j)(2,1,2,0)");
  test_for_zero(t4ddg_1(2,2,1,1) - (t2_1(2,1))
		,"T4ddg(Num,i,Num,j)(2,1,2,1)");
  test_for_zero(t4ddg_1(2,2,1,2) - (t2_1(2,2))
		,"T4ddg(Num,i,Num,j)(2,1,2,2)");

  t2_1(i,j)=t4ddg_1(2,i,2,j);
  test_for_zero(t4ddg_1(2,0,2,0) - (t2_1(0,0))
		,"T4ddg(Num,i,Num,j)(2,2,0,0)");
  test_for_zero(t4ddg_1(2,0,2,1) - (t2_1(0,1))
		,"T4ddg(Num,i,Num,j)(2,2,0,1)");
  test_for_zero(t4ddg_1(2,0,2,2) - (t2_1(0,2))
		,"T4ddg(Num,i,Num,j)(2,2,0,2)");
  test_for_zero(t4ddg_1(2,1,2,0) - (t2_1(1,0))
		,"T4ddg(Num,i,Num,j)(2,2,1,0)");
  test_for_zero(t4ddg_1(2,1,2,1) - (t2_1(1,1))
		,"T4ddg(Num,i,Num,j)(2,2,1,1)");
  test_for_zero(t4ddg_1(2,1,2,2) - (t2_1(1,2))
		,"T4ddg(Num,i,Num,j)(2,2,1,2)");
  test_for_zero(t4ddg_1(2,2,2,0) - (t2_1(2,0))
		,"T4ddg(Num,i,Num,j)(2,2,2,0)");
  test_for_zero(t4ddg_1(2,2,2,1) - (t2_1(2,1))
		,"T4ddg(Num,i,Num,j)(2,2,2,1)");
  test_for_zero(t4ddg_1(2,2,2,2) - (t2_1(2,2))
		,"T4ddg(Num,i,Num,j)(2,2,2,2)");

  t3dg_1(j,k,i)=t4ddg_1(0,i,j,k);
  test_for_zero(t3dg_1(0,0,0) - t4ddg_1(0,0,0,0)
		,"T4ddg(Num,i,j,k)(0,0,0,0)");
  test_for_zero(t3dg_1(0,1,0) - t4ddg_1(0,0,0,1)
		,"T4ddg(Num,i,j,k)(0,0,0,1)");
  test_for_zero(t3dg_1(0,2,0) - t4ddg_1(0,0,0,2)
		,"T4ddg(Num,i,j,k)(0,0,0,2)");
  test_for_zero(t3dg_1(1,0,0) - t4ddg_1(0,0,1,0)
		,"T4ddg(Num,i,j,k)(0,0,1,0)");
  test_for_zero(t3dg_1(1,1,0) - t4ddg_1(0,0,1,1)
		,"T4ddg(Num,i,j,k)(0,0,1,1)");
  test_for_zero(t3dg_1(1,2,0) - t4ddg_1(0,0,1,2)
		,"T4ddg(Num,i,j,k)(0,0,1,2)");
  test_for_zero(t3dg_1(2,0,0) - t4ddg_1(0,0,2,0)
		,"T4ddg(Num,i,j,k)(0,0,2,0)");
  test_for_zero(t3dg_1(2,1,0) - t4ddg_1(0,0,2,1)
		,"T4ddg(Num,i,j,k)(0,0,2,1)");
  test_for_zero(t3dg_1(2,2,0) - t4ddg_1(0,0,2,2)
		,"T4ddg(Num,i,j,k)(0,0,2,2)");
  test_for_zero(t3dg_1(0,0,1) - t4ddg_1(0,1,0,0)
		,"T4ddg(Num,i,j,k)(0,1,0,0)");
  test_for_zero(t3dg_1(0,1,1) - t4ddg_1(0,1,0,1)
		,"T4ddg(Num,i,j,k)(0,1,0,1)");
  test_for_zero(t3dg_1(0,2,1) - t4ddg_1(0,1,0,2)
		,"T4ddg(Num,i,j,k)(0,1,0,2)");
  test_for_zero(t3dg_1(1,0,1) - t4ddg_1(0,1,1,0)
		,"T4ddg(Num,i,j,k)(0,1,1,0)");
  test_for_zero(t3dg_1(1,1,1) - t4ddg_1(0,1,1,1)
		,"T4ddg(Num,i,j,k)(0,1,1,1)");
  test_for_zero(t3dg_1(1,2,1) - t4ddg_1(0,1,1,2)
		,"T4ddg(Num,i,j,k)(0,1,1,2)");
  test_for_zero(t3dg_1(2,0,1) - t4ddg_1(0,1,2,0)
		,"T4ddg(Num,i,j,k)(0,1,2,0)");
  test_for_zero(t3dg_1(2,1,1) - t4ddg_1(0,1,2,1)
		,"T4ddg(Num,i,j,k)(0,1,2,1)");
  test_for_zero(t3dg_1(2,2,1) - t4ddg_1(0,1,2,2)
		,"T4ddg(Num,i,j,k)(0,1,2,2)");
  test_for_zero(t3dg_1(0,0,2) - t4ddg_1(0,2,0,0)
		,"T4ddg(Num,i,j,k)(0,2,0,0)");
  test_for_zero(t3dg_1(0,1,2) - t4ddg_1(0,2,0,1)
		,"T4ddg(Num,i,j,k)(0,2,0,1)");
  test_for_zero(t3dg_1(0,2,2) - t4ddg_1(0,2,0,2)
		,"T4ddg(Num,i,j,k)(0,2,0,2)");
  test_for_zero(t3dg_1(1,0,2) - t4ddg_1(0,2,1,0)
		,"T4ddg(Num,i,j,k)(0,2,1,0)");
  test_for_zero(t3dg_1(1,1,2) - t4ddg_1(0,2,1,1)
		,"T4ddg(Num,i,j,k)(0,2,1,1)");
  test_for_zero(t3dg_1(1,2,2) - t4ddg_1(0,2,1,2)
		,"T4ddg(Num,i,j,k)(0,2,1,2)");
  test_for_zero(t3dg_1(2,0,2) - t4ddg_1(0,2,2,0)
		,"T4ddg(Num,i,j,k)(0,2,2,0)");
  test_for_zero(t3dg_1(2,1,2) - t4ddg_1(0,2,2,1)
		,"T4ddg(Num,i,j,k)(0,2,2,1)");
  test_for_zero(t3dg_1(2,2,2) - t4ddg_1(0,2,2,2)
		,"T4ddg(Num,i,j,k)(0,2,2,2)");

  t3dg_1(j,k,i)=t4ddg_1(1,i,j,k);
  test_for_zero(t3dg_1(0,0,0) - t4ddg_1(1,0,0,0)
		,"T4ddg(Num,i,j,k)(1,0,0,0)");
  test_for_zero(t3dg_1(0,1,0) - t4ddg_1(1,0,0,1)
		,"T4ddg(Num,i,j,k)(1,0,0,1)");
  test_for_zero(t3dg_1(0,2,0) - t4ddg_1(1,0,0,2)
		,"T4ddg(Num,i,j,k)(1,0,0,2)");
  test_for_zero(t3dg_1(1,0,0) - t4ddg_1(1,0,1,0)
		,"T4ddg(Num,i,j,k)(1,0,1,0)");
  test_for_zero(t3dg_1(1,1,0) - t4ddg_1(1,0,1,1)
		,"T4ddg(Num,i,j,k)(1,0,1,1)");
  test_for_zero(t3dg_1(1,2,0) - t4ddg_1(1,0,1,2)
		,"T4ddg(Num,i,j,k)(1,0,1,2)");
  test_for_zero(t3dg_1(2,0,0) - t4ddg_1(1,0,2,0)
		,"T4ddg(Num,i,j,k)(1,0,2,0)");
  test_for_zero(t3dg_1(2,1,0) - t4ddg_1(1,0,2,1)
		,"T4ddg(Num,i,j,k)(1,0,2,1)");
  test_for_zero(t3dg_1(2,2,0) - t4ddg_1(1,0,2,2)
		,"T4ddg(Num,i,j,k)(1,0,2,2)");
  test_for_zero(t3dg_1(0,0,1) - t4ddg_1(1,1,0,0)
		,"T4ddg(Num,i,j,k)(1,1,0,0)");
  test_for_zero(t3dg_1(0,1,1) - t4ddg_1(1,1,0,1)
		,"T4ddg(Num,i,j,k)(1,1,0,1)");
  test_for_zero(t3dg_1(0,2,1) - t4ddg_1(1,1,0,2)
		,"T4ddg(Num,i,j,k)(1,1,0,2)");
  test_for_zero(t3dg_1(1,0,1) - t4ddg_1(1,1,1,0)
		,"T4ddg(Num,i,j,k)(1,1,1,0)");
  test_for_zero(t3dg_1(1,1,1) - t4ddg_1(1,1,1,1)
		,"T4ddg(Num,i,j,k)(1,1,1,1)");
  test_for_zero(t3dg_1(1,2,1) - t4ddg_1(1,1,1,2)
		,"T4ddg(Num,i,j,k)(1,1,1,2)");
  test_for_zero(t3dg_1(2,0,1) - t4ddg_1(1,1,2,0)
		,"T4ddg(Num,i,j,k)(1,1,2,0)");
  test_for_zero(t3dg_1(2,1,1) - t4ddg_1(1,1,2,1)
		,"T4ddg(Num,i,j,k)(1,1,2,1)");
  test_for_zero(t3dg_1(2,2,1) - t4ddg_1(1,1,2,2)
		,"T4ddg(Num,i,j,k)(1,1,2,2)");
  test_for_zero(t3dg_1(0,0,2) - t4ddg_1(1,2,0,0)
		,"T4ddg(Num,i,j,k)(1,2,0,0)");
  test_for_zero(t3dg_1(0,1,2) - t4ddg_1(1,2,0,1)
		,"T4ddg(Num,i,j,k)(1,2,0,1)");
  test_for_zero(t3dg_1(0,2,2) - t4ddg_1(1,2,0,2)
		,"T4ddg(Num,i,j,k)(1,2,0,2)");
  test_for_zero(t3dg_1(1,0,2) - t4ddg_1(1,2,1,0)
		,"T4ddg(Num,i,j,k)(1,2,1,0)");
  test_for_zero(t3dg_1(1,1,2) - t4ddg_1(1,2,1,1)
		,"T4ddg(Num,i,j,k)(1,2,1,1)");
  test_for_zero(t3dg_1(1,2,2) - t4ddg_1(1,2,1,2)
		,"T4ddg(Num,i,j,k)(1,2,1,2)");
  test_for_zero(t3dg_1(2,0,2) - t4ddg_1(1,2,2,0)
		,"T4ddg(Num,i,j,k)(1,2,2,0)");
  test_for_zero(t3dg_1(2,1,2) - t4ddg_1(1,2,2,1)
		,"T4ddg(Num,i,j,k)(1,2,2,1)");
  test_for_zero(t3dg_1(2,2,2) - t4ddg_1(1,2,2,2)
		,"T4ddg(Num,i,j,k)(1,2,2,2)");

  t3dg_1(j,k,i)=t4ddg_1(2,i,j,k);
  test_for_zero(t3dg_1(0,0,0) - t4ddg_1(2,0,0,0)
		,"T4ddg(Num,i,j,k)(2,0,0,0)");
  test_for_zero(t3dg_1(0,1,0) - t4ddg_1(2,0,0,1)
		,"T4ddg(Num,i,j,k)(2,0,0,1)");
  test_for_zero(t3dg_1(0,2,0) - t4ddg_1(2,0,0,2)
		,"T4ddg(Num,i,j,k)(2,0,0,2)");
  test_for_zero(t3dg_1(1,0,0) - t4ddg_1(2,0,1,0)
		,"T4ddg(Num,i,j,k)(2,0,1,0)");
  test_for_zero(t3dg_1(1,1,0) - t4ddg_1(2,0,1,1)
		,"T4ddg(Num,i,j,k)(2,0,1,1)");
  test_for_zero(t3dg_1(1,2,0) - t4ddg_1(2,0,1,2)
		,"T4ddg(Num,i,j,k)(2,0,1,2)");
  test_for_zero(t3dg_1(2,0,0) - t4ddg_1(2,0,2,0)
		,"T4ddg(Num,i,j,k)(2,0,2,0)");
  test_for_zero(t3dg_1(2,1,0) - t4ddg_1(2,0,2,1)
		,"T4ddg(Num,i,j,k)(2,0,2,1)");
  test_for_zero(t3dg_1(2,2,0) - t4ddg_1(2,0,2,2)
		,"T4ddg(Num,i,j,k)(2,0,2,2)");
  test_for_zero(t3dg_1(0,0,1) - t4ddg_1(2,1,0,0)
		,"T4ddg(Num,i,j,k)(2,1,0,0)");
  test_for_zero(t3dg_1(0,1,1) - t4ddg_1(2,1,0,1)
		,"T4ddg(Num,i,j,k)(2,1,0,1)");
  test_for_zero(t3dg_1(0,2,1) - t4ddg_1(2,1,0,2)
		,"T4ddg(Num,i,j,k)(2,1,0,2)");
  test_for_zero(t3dg_1(1,0,1) - t4ddg_1(2,1,1,0)
		,"T4ddg(Num,i,j,k)(2,1,1,0)");
  test_for_zero(t3dg_1(1,1,1) - t4ddg_1(2,1,1,1)
		,"T4ddg(Num,i,j,k)(2,1,1,1)");
  test_for_zero(t3dg_1(1,2,1) - t4ddg_1(2,1,1,2)
		,"T4ddg(Num,i,j,k)(2,1,1,2)");
  test_for_zero(t3dg_1(2,0,1) - t4ddg_1(2,1,2,0)
		,"T4ddg(Num,i,j,k)(2,1,2,0)");
  test_for_zero(t3dg_1(2,1,1) - t4ddg_1(2,1,2,1)
		,"T4ddg(Num,i,j,k)(2,1,2,1)");
  test_for_zero(t3dg_1(2,2,1) - t4ddg_1(2,1,2,2)
		,"T4ddg(Num,i,j,k)(2,1,2,2)");
  test_for_zero(t3dg_1(0,0,2) - t4ddg_1(2,2,0,0)
		,"T4ddg(Num,i,j,k)(2,2,0,0)");
  test_for_zero(t3dg_1(0,1,2) - t4ddg_1(2,2,0,1)
		,"T4ddg(Num,i,j,k)(2,2,0,1)");
  test_for_zero(t3dg_1(0,2,2) - t4ddg_1(2,2,0,2)
		,"T4ddg(Num,i,j,k)(2,2,0,2)");
  test_for_zero(t3dg_1(1,0,2) - t4ddg_1(2,2,1,0)
		,"T4ddg(Num,i,j,k)(2,2,1,0)");
  test_for_zero(t3dg_1(1,1,2) - t4ddg_1(2,2,1,1)
		,"T4ddg(Num,i,j,k)(2,2,1,1)");
  test_for_zero(t3dg_1(1,2,2) - t4ddg_1(2,2,1,2)
		,"T4ddg(Num,i,j,k)(2,2,1,2)");
  test_for_zero(t3dg_1(2,0,2) - t4ddg_1(2,2,2,0)
		,"T4ddg(Num,i,j,k)(2,2,2,0)");
  test_for_zero(t3dg_1(2,1,2) - t4ddg_1(2,2,2,1)
		,"T4ddg(Num,i,j,k)(2,2,2,1)");
  test_for_zero(t3dg_1(2,2,2) - t4ddg_1(2,2,2,2)
		,"T4ddg(Num,i,j,k)(2,2,2,2)");

  /* Only one index int the first slot*/

  t1_1(i)=t4ddg_1(i,0,0,0);
  test_for_zero(t1_1(0) - t4ddg_1(0,0,0,0)
		,"T4ddg(i,Num,Num,Num)(0,0,0,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,0,0)
		,"T4ddg(i,Num,Num,Num)(0,0,0,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,0,0)
		,"T4ddg(i,Num,Num,Num)(0,0,0,2)");
  t1_1(i)=t4ddg_1(i,0,0,1);
  test_for_zero(t1_1(0) - t4ddg_1(0,0,0,1)
		,"T4ddg(i,Num,Num,Num)(0,0,1,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,0,1)
		,"T4ddg(i,Num,Num,Num)(0,0,1,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,0,1)
		,"T4ddg(i,Num,Num,Num)(0,0,1,2)");
  t1_1(i)=t4ddg_1(i,0,0,2);
  test_for_zero(t1_1(0) - t4ddg_1(0,0,0,2)
		,"T4ddg(i,Num,Num,Num)(0,0,2,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,0,2)
		,"T4ddg(i,Num,Num,Num)(0,0,2,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,0,2)
		,"T4ddg(i,Num,Num,Num)(0,0,2,2)");
  t1_1(i)=t4ddg_1(i,0,1,0);
  test_for_zero(t1_1(0) - t4ddg_1(0,0,1,0)
		,"T4ddg(i,Num,Num,Num)(0,1,0,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,1,0)
		,"T4ddg(i,Num,Num,Num)(0,1,0,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,1,0)
		,"T4ddg(i,Num,Num,Num)(0,1,0,2)");
  t1_1(i)=t4ddg_1(i,0,1,1);
  test_for_zero(t1_1(0) - t4ddg_1(0,0,1,1)
		,"T4ddg(i,Num,Num,Num)(0,1,1,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,1,1)
		,"T4ddg(i,Num,Num,Num)(0,1,1,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,1,1)
		,"T4ddg(i,Num,Num,Num)(0,1,1,2)");
  t1_1(i)=t4ddg_1(i,0,1,2);
  test_for_zero(t1_1(0) - t4ddg_1(0,0,1,2)
		,"T4ddg(i,Num,Num,Num)(0,1,2,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,1,2)
		,"T4ddg(i,Num,Num,Num)(0,1,2,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,1,2)
		,"T4ddg(i,Num,Num,Num)(0,1,2,2)");
  t1_1(i)=t4ddg_1(i,0,2,0);
  test_for_zero(t1_1(0) - t4ddg_1(0,0,2,0)
		,"T4ddg(i,Num,Num,Num)(0,2,0,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,2,0)
		,"T4ddg(i,Num,Num,Num)(0,2,0,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,2,0)
		,"T4ddg(i,Num,Num,Num)(0,2,0,2)");
  t1_1(i)=t4ddg_1(i,0,2,1);
  test_for_zero(t1_1(0) - t4ddg_1(0,0,2,1)
		,"T4ddg(i,Num,Num,Num)(0,2,1,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,2,1)
		,"T4ddg(i,Num,Num,Num)(0,2,1,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,2,1)
		,"T4ddg(i,Num,Num,Num)(0,2,1,2)");
  t1_1(i)=t4ddg_1(i,0,2,2);
  test_for_zero(t1_1(0) - t4ddg_1(0,0,2,2)
		,"T4ddg(i,Num,Num,Num)(0,2,2,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,2,2)
		,"T4ddg(i,Num,Num,Num)(0,2,2,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,2,2)
		,"T4ddg(i,Num,Num,Num)(0,2,2,2)");
  t1_1(i)=t4ddg_1(i,1,0,0);
  test_for_zero(t1_1(0) - t4ddg_1(0,1,0,0)
		,"T4ddg(i,Num,Num,Num)(1,0,0,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,1,0,0)
		,"T4ddg(i,Num,Num,Num)(1,0,0,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,1,0,0)
		,"T4ddg(i,Num,Num,Num)(1,0,0,2)");
  t1_1(i)=t4ddg_1(i,1,0,1);
  test_for_zero(t1_1(0) - t4ddg_1(0,1,0,1)
		,"T4ddg(i,Num,Num,Num)(1,0,1,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,1,0,1)
		,"T4ddg(i,Num,Num,Num)(1,0,1,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,1,0,1)
		,"T4ddg(i,Num,Num,Num)(1,0,1,2)");
  t1_1(i)=t4ddg_1(i,1,0,2);
  test_for_zero(t1_1(0) - t4ddg_1(0,1,0,2)
		,"T4ddg(i,Num,Num,Num)(1,0,2,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,1,0,2)
		,"T4ddg(i,Num,Num,Num)(1,0,2,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,1,0,2)
		,"T4ddg(i,Num,Num,Num)(1,0,2,2)");
  t1_1(i)=t4ddg_1(i,1,1,0);
  test_for_zero(t1_1(0) - t4ddg_1(0,1,1,0)
		,"T4ddg(i,Num,Num,Num)(1,1,0,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,1,1,0)
		,"T4ddg(i,Num,Num,Num)(1,1,0,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,1,1,0)
		,"T4ddg(i,Num,Num,Num)(1,1,0,2)");
  t1_1(i)=t4ddg_1(i,1,1,1);
  test_for_zero(t1_1(0) - t4ddg_1(0,1,1,1)
		,"T4ddg(i,Num,Num,Num)(1,1,1,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,1,1,1)
		,"T4ddg(i,Num,Num,Num)(1,1,1,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,1,1,1)
		,"T4ddg(i,Num,Num,Num)(1,1,1,2)");
  t1_1(i)=t4ddg_1(i,1,1,2);
  test_for_zero(t1_1(0) - t4ddg_1(0,1,1,2)
		,"T4ddg(i,Num,Num,Num)(1,1,2,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,1,1,2)
		,"T4ddg(i,Num,Num,Num)(1,1,2,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,1,1,2)
		,"T4ddg(i,Num,Num,Num)(1,1,2,2)");
  t1_1(i)=t4ddg_1(i,1,2,0);
  test_for_zero(t1_1(0) - t4ddg_1(0,1,2,0)
		,"T4ddg(i,Num,Num,Num)(1,2,0,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,1,2,0)
		,"T4ddg(i,Num,Num,Num)(1,2,0,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,1,2,0)
		,"T4ddg(i,Num,Num,Num)(1,2,0,2)");
  t1_1(i)=t4ddg_1(i,1,2,1);
  test_for_zero(t1_1(0) - t4ddg_1(0,1,2,1)
		,"T4ddg(i,Num,Num,Num)(1,2,1,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,1,2,1)
		,"T4ddg(i,Num,Num,Num)(1,2,1,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,1,2,1)
		,"T4ddg(i,Num,Num,Num)(1,2,1,2)");
  t1_1(i)=t4ddg_1(i,1,2,2);
  test_for_zero(t1_1(0) - t4ddg_1(0,1,2,2)
		,"T4ddg(i,Num,Num,Num)(1,2,2,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,1,2,2)
		,"T4ddg(i,Num,Num,Num)(1,2,2,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,1,2,2)
		,"T4ddg(i,Num,Num,Num)(1,2,2,2)");
  t1_1(i)=t4ddg_1(i,2,0,0);
  test_for_zero(t1_1(0) - t4ddg_1(0,2,0,0)
		,"T4ddg(i,Num,Num,Num)(2,0,0,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,2,0,0)
		,"T4ddg(i,Num,Num,Num)(2,0,0,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,2,0,0)
		,"T4ddg(i,Num,Num,Num)(2,0,0,2)");
  t1_1(i)=t4ddg_1(i,2,0,1);
  test_for_zero(t1_1(0) - t4ddg_1(0,2,0,1)
		,"T4ddg(i,Num,Num,Num)(2,0,1,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,2,0,1)
		,"T4ddg(i,Num,Num,Num)(2,0,1,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,2,0,1)
		,"T4ddg(i,Num,Num,Num)(2,0,1,2)");
  t1_1(i)=t4ddg_1(i,2,0,2);
  test_for_zero(t1_1(0) - t4ddg_1(0,2,0,2)
		,"T4ddg(i,Num,Num,Num)(2,0,2,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,2,0,2)
		,"T4ddg(i,Num,Num,Num)(2,0,2,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,2,0,2)
		,"T4ddg(i,Num,Num,Num)(2,0,2,2)");
  t1_1(i)=t4ddg_1(i,2,1,0);
  test_for_zero(t1_1(0) - t4ddg_1(0,2,1,0)
		,"T4ddg(i,Num,Num,Num)(2,1,0,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,2,1,0)
		,"T4ddg(i,Num,Num,Num)(2,1,0,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,2,1,0)
		,"T4ddg(i,Num,Num,Num)(2,1,0,2)");
  t1_1(i)=t4ddg_1(i,2,1,1);
  test_for_zero(t1_1(0) - t4ddg_1(0,2,1,1)
		,"T4ddg(i,Num,Num,Num)(2,1,1,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,2,1,1)
		,"T4ddg(i,Num,Num,Num)(2,1,1,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,2,1,1)
		,"T4ddg(i,Num,Num,Num)(2,1,1,2)");
  t1_1(i)=t4ddg_1(i,2,1,2);
  test_for_zero(t1_1(0) - t4ddg_1(0,2,1,2)
		,"T4ddg(i,Num,Num,Num)(2,1,2,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,2,1,2)
		,"T4ddg(i,Num,Num,Num)(2,1,2,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,2,1,2)
		,"T4ddg(i,Num,Num,Num)(2,1,2,2)");
  t1_1(i)=t4ddg_1(i,2,2,0);
  test_for_zero(t1_1(0) - t4ddg_1(0,2,2,0)
		,"T4ddg(i,Num,Num,Num)(2,2,0,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,2,2,0)
		,"T4ddg(i,Num,Num,Num)(2,2,0,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,2,2,0)
		,"T4ddg(i,Num,Num,Num)(2,2,0,2)");
  t1_1(i)=t4ddg_1(i,2,2,1);
  test_for_zero(t1_1(0) - t4ddg_1(0,2,2,1)
		,"T4ddg(i,Num,Num,Num)(2,2,1,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,2,2,1)
		,"T4ddg(i,Num,Num,Num)(2,2,1,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,2,2,1)
		,"T4ddg(i,Num,Num,Num)(2,2,1,2)");
  t1_1(i)=t4ddg_1(i,2,2,2);
  test_for_zero(t1_1(0) - t4ddg_1(0,2,2,2)
		,"T4ddg(i,Num,Num,Num)(2,2,2,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,2,2,2)
		,"T4ddg(i,Num,Num,Num)(2,2,2,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,2,2,2)
		,"T4ddg(i,Num,Num,Num)(2,2,2,2)");

  /* Only one index in the second slot */

  t1_1(i)=t4ddg_1(0,i,0,0);
  test_for_zero(t1_1(0) - t4ddg_1(0,0,0,0)
		,"T4ddg(Num,i,Num,Num)(0,0,0,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,0,0)
		,"T4ddg(Num,i,Num,Num)(0,0,0,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,0,0)
		,"T4ddg(Num,i,Num,Num)(0,0,0,2)");
  t1_1(i)=t4ddg_1(0,i,0,1);
  test_for_zero(t1_1(0) - t4ddg_1(0,0,0,1)
		,"T4ddg(Num,i,Num,Num)(0,0,1,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,0,1)
		,"T4ddg(Num,i,Num,Num)(0,0,1,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,0,1)
		,"T4ddg(Num,i,Num,Num)(0,0,1,2)");
  t1_1(i)=t4ddg_1(0,i,0,2);
  test_for_zero(t1_1(0) - t4ddg_1(0,0,0,2)
		,"T4ddg(Num,i,Num,Num)(0,0,2,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,0,2)
		,"T4ddg(Num,i,Num,Num)(0,0,2,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,0,2)
		,"T4ddg(Num,i,Num,Num)(0,0,2,2)");
  t1_1(i)=t4ddg_1(0,i,1,0);
  test_for_zero(t1_1(0) - t4ddg_1(0,0,1,0)
		,"T4ddg(Num,i,Num,Num)(0,1,0,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,1,0)
		,"T4ddg(Num,i,Num,Num)(0,1,0,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,1,0)
		,"T4ddg(Num,i,Num,Num)(0,1,0,2)");
  t1_1(i)=t4ddg_1(0,i,1,1);
  test_for_zero(t1_1(0) - t4ddg_1(0,0,1,1)
		,"T4ddg(Num,i,Num,Num)(0,1,1,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,1,1)
		,"T4ddg(Num,i,Num,Num)(0,1,1,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,1,1)
		,"T4ddg(Num,i,Num,Num)(0,1,1,2)");
  t1_1(i)=t4ddg_1(0,i,1,2);
  test_for_zero(t1_1(0) - t4ddg_1(0,0,1,2)
		,"T4ddg(Num,i,Num,Num)(0,1,2,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,1,2)
		,"T4ddg(Num,i,Num,Num)(0,1,2,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,1,2)
		,"T4ddg(Num,i,Num,Num)(0,1,2,2)");
  t1_1(i)=t4ddg_1(0,i,2,0);
  test_for_zero(t1_1(0) - t4ddg_1(0,0,2,0)
		,"T4ddg(Num,i,Num,Num)(0,2,0,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,2,0)
		,"T4ddg(Num,i,Num,Num)(0,2,0,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,2,0)
		,"T4ddg(Num,i,Num,Num)(0,2,0,2)");
  t1_1(i)=t4ddg_1(0,i,2,1);
  test_for_zero(t1_1(0) - t4ddg_1(0,0,2,1)
		,"T4ddg(Num,i,Num,Num)(0,2,1,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,2,1)
		,"T4ddg(Num,i,Num,Num)(0,2,1,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,2,1)
		,"T4ddg(Num,i,Num,Num)(0,2,1,2)");
  t1_1(i)=t4ddg_1(0,i,2,2);
  test_for_zero(t1_1(0) - t4ddg_1(0,0,2,2)
		,"T4ddg(Num,i,Num,Num)(0,2,2,0)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,2,2)
		,"T4ddg(Num,i,Num,Num)(0,2,2,1)");
  test_for_zero(t1_1(0) - t4ddg_1(0,0,2,2)
		,"T4ddg(Num,i,Num,Num)(0,2,2,2)");
  t1_1(i)=t4ddg_1(1,i,0,0);
  test_for_zero(t1_1(0) - t4ddg_1(1,0,0,0)
		,"T4ddg(Num,i,Num,Num)(1,0,0,0)");
  test_for_zero(t1_1(0) - t4ddg_1(1,0,0,0)
		,"T4ddg(Num,i,Num,Num)(1,0,0,1)");
  test_for_zero(t1_1(0) - t4ddg_1(1,0,0,0)
		,"T4ddg(Num,i,Num,Num)(1,0,0,2)");
  t1_1(i)=t4ddg_1(1,i,0,1);
  test_for_zero(t1_1(0) - t4ddg_1(1,0,0,1)
		,"T4ddg(Num,i,Num,Num)(1,0,1,0)");
  test_for_zero(t1_1(0) - t4ddg_1(1,0,0,1)
		,"T4ddg(Num,i,Num,Num)(1,0,1,1)");
  test_for_zero(t1_1(0) - t4ddg_1(1,0,0,1)
		,"T4ddg(Num,i,Num,Num)(1,0,1,2)");
  t1_1(i)=t4ddg_1(1,i,0,2);
  test_for_zero(t1_1(0) - t4ddg_1(1,0,0,2)
		,"T4ddg(Num,i,Num,Num)(1,0,2,0)");
  test_for_zero(t1_1(0) - t4ddg_1(1,0,0,2)
		,"T4ddg(Num,i,Num,Num)(1,0,2,1)");
  test_for_zero(t1_1(0) - t4ddg_1(1,0,0,2)
		,"T4ddg(Num,i,Num,Num)(1,0,2,2)");
  t1_1(i)=t4ddg_1(1,i,1,0);
  test_for_zero(t1_1(0) - t4ddg_1(1,0,1,0)
		,"T4ddg(Num,i,Num,Num)(1,1,0,0)");
  test_for_zero(t1_1(0) - t4ddg_1(1,0,1,0)
		,"T4ddg(Num,i,Num,Num)(1,1,0,1)");
  test_for_zero(t1_1(0) - t4ddg_1(1,0,1,0)
		,"T4ddg(Num,i,Num,Num)(1,1,0,2)");
  t1_1(i)=t4ddg_1(1,i,1,1);
  test_for_zero(t1_1(0) - t4ddg_1(1,0,1,1)
		,"T4ddg(Num,i,Num,Num)(1,1,1,0)");
  test_for_zero(t1_1(0) - t4ddg_1(1,0,1,1)
		,"T4ddg(Num,i,Num,Num)(1,1,1,1)");
  test_for_zero(t1_1(0) - t4ddg_1(1,0,1,1)
		,"T4ddg(Num,i,Num,Num)(1,1,1,2)");
  t1_1(i)=t4ddg_1(1,i,1,2);
  test_for_zero(t1_1(0) - t4ddg_1(1,0,1,2)
		,"T4ddg(Num,i,Num,Num)(1,1,2,0)");
  test_for_zero(t1_1(0) - t4ddg_1(1,0,1,2)
		,"T4ddg(Num,i,Num,Num)(1,1,2,1)");
  test_for_zero(t1_1(0) - t4ddg_1(1,0,1,2)
		,"T4ddg(Num,i,Num,Num)(1,1,2,2)");
  t1_1(i)=t4ddg_1(1,i,2,0);
  test_for_zero(t1_1(0) - t4ddg_1(1,0,2,0)
		,"T4ddg(Num,i,Num,Num)(1,2,0,0)");
  test_for_zero(t1_1(0) - t4ddg_1(1,0,2,0)
		,"T4ddg(Num,i,Num,Num)(1,2,0,1)");
  test_for_zero(t1_1(0) - t4ddg_1(1,0,2,0)
		,"T4ddg(Num,i,Num,Num)(1,2,0,2)");
  t1_1(i)=t4ddg_1(1,i,2,1);
  test_for_zero(t1_1(0) - t4ddg_1(1,0,2,1)
		,"T4ddg(Num,i,Num,Num)(1,2,1,0)");
  test_for_zero(t1_1(0) - t4ddg_1(1,0,2,1)
		,"T4ddg(Num,i,Num,Num)(1,2,1,1)");
  test_for_zero(t1_1(0) - t4ddg_1(1,0,2,1)
		,"T4ddg(Num,i,Num,Num)(1,2,1,2)");
  t1_1(i)=t4ddg_1(1,i,2,2);
  test_for_zero(t1_1(0) - t4ddg_1(1,0,2,2)
		,"T4ddg(Num,i,Num,Num)(1,2,2,0)");
  test_for_zero(t1_1(0) - t4ddg_1(1,0,2,2)
		,"T4ddg(Num,i,Num,Num)(1,2,2,1)");
  test_for_zero(t1_1(0) - t4ddg_1(1,0,2,2)
		,"T4ddg(Num,i,Num,Num)(1,2,2,2)");
  t1_1(i)=t4ddg_1(2,i,0,0);
  test_for_zero(t1_1(0) - t4ddg_1(2,0,0,0)
		,"T4ddg(Num,i,Num,Num)(2,0,0,0)");
  test_for_zero(t1_1(0) - t4ddg_1(2,0,0,0)
		,"T4ddg(Num,i,Num,Num)(2,0,0,1)");
  test_for_zero(t1_1(0) - t4ddg_1(2,0,0,0)
		,"T4ddg(Num,i,Num,Num)(2,0,0,2)");
  t1_1(i)=t4ddg_1(2,i,0,1);
  test_for_zero(t1_1(0) - t4ddg_1(2,0,0,1)
		,"T4ddg(Num,i,Num,Num)(2,0,1,0)");
  test_for_zero(t1_1(0) - t4ddg_1(2,0,0,1)
		,"T4ddg(Num,i,Num,Num)(2,0,1,1)");
  test_for_zero(t1_1(0) - t4ddg_1(2,0,0,1)
		,"T4ddg(Num,i,Num,Num)(2,0,1,2)");
  t1_1(i)=t4ddg_1(2,i,0,2);
  test_for_zero(t1_1(0) - t4ddg_1(2,0,0,2)
		,"T4ddg(Num,i,Num,Num)(2,0,2,0)");
  test_for_zero(t1_1(0) - t4ddg_1(2,0,0,2)
		,"T4ddg(Num,i,Num,Num)(2,0,2,1)");
  test_for_zero(t1_1(0) - t4ddg_1(2,0,0,2)
		,"T4ddg(Num,i,Num,Num)(2,0,2,2)");
  t1_1(i)=t4ddg_1(2,i,1,0);
  test_for_zero(t1_1(0) - t4ddg_1(2,0,1,0)
		,"T4ddg(Num,i,Num,Num)(2,1,0,0)");
  test_for_zero(t1_1(0) - t4ddg_1(2,0,1,0)
		,"T4ddg(Num,i,Num,Num)(2,1,0,1)");
  test_for_zero(t1_1(0) - t4ddg_1(2,0,1,0)
		,"T4ddg(Num,i,Num,Num)(2,1,0,2)");
  t1_1(i)=t4ddg_1(2,i,1,1);
  test_for_zero(t1_1(0) - t4ddg_1(2,0,1,1)
		,"T4ddg(Num,i,Num,Num)(2,1,1,0)");
  test_for_zero(t1_1(0) - t4ddg_1(2,0,1,1)
		,"T4ddg(Num,i,Num,Num)(2,1,1,1)");
  test_for_zero(t1_1(0) - t4ddg_1(2,0,1,1)
		,"T4ddg(Num,i,Num,Num)(2,1,1,2)");
  t1_1(i)=t4ddg_1(2,i,1,2);
  test_for_zero(t1_1(0) - t4ddg_1(2,0,1,2)
		,"T4ddg(Num,i,Num,Num)(2,1,2,0)");
  test_for_zero(t1_1(0) - t4ddg_1(2,0,1,2)
		,"T4ddg(Num,i,Num,Num)(2,1,2,1)");
  test_for_zero(t1_1(0) - t4ddg_1(2,0,1,2)
		,"T4ddg(Num,i,Num,Num)(2,1,2,2)");
  t1_1(i)=t4ddg_1(2,i,2,0);
  test_for_zero(t1_1(0) - t4ddg_1(2,0,2,0)
		,"T4ddg(Num,i,Num,Num)(2,2,0,0)");
  test_for_zero(t1_1(0) - t4ddg_1(2,0,2,0)
		,"T4ddg(Num,i,Num,Num)(2,2,0,1)");
  test_for_zero(t1_1(0) - t4ddg_1(2,0,2,0)
		,"T4ddg(Num,i,Num,Num)(2,2,0,2)");
  t1_1(i)=t4ddg_1(2,i,2,1);
  test_for_zero(t1_1(0) - t4ddg_1(2,0,2,1)
		,"T4ddg(Num,i,Num,Num)(2,2,1,0)");
  test_for_zero(t1_1(0) - t4ddg_1(2,0,2,1)
		,"T4ddg(Num,i,Num,Num)(2,2,1,1)");
  test_for_zero(t1_1(0) - t4ddg_1(2,0,2,1)
		,"T4ddg(Num,i,Num,Num)(2,2,1,2)");
  t1_1(i)=t4ddg_1(2,i,2,2);
  test_for_zero(t1_1(0) - t4ddg_1(2,0,2,2)
		,"T4ddg(Num,i,Num,Num)(2,2,2,0)");
  test_for_zero(t1_1(0) - t4ddg_1(2,0,2,2)
		,"T4ddg(Num,i,Num,Num)(2,2,2,1)");
  test_for_zero(t1_1(0) - t4ddg_1(2,0,2,2)
		,"T4ddg(Num,i,Num,Num)(2,2,2,2)");

  cout << endl;
}
