#include <iostream>
#include "../../../FTensor.h"
#include "../test_for_zero.h"
using namespace FTensor;
using namespace std;

void test_T3asIX(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	       Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	       Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	       const Tensor2<double,3,3> &t2_3,
	       Tensor2_symmetric<double,3> &t2s_1,
	       const Tensor2_symmetric<double,3> &t2s_2,
	       const Tensor2_symmetric<double,3> &t2s_3,
	       Tensor3_dg<double,3,3> &t3dg_1,
	       const Tensor3_dg<double,3,3> &t3dg_2,
	       const Tensor3_dg<double,3,3> &t3dg_3,
	       Tensor3_antisymmetric<double,3,3> &t3as_1,
	       const Tensor3_antisymmetric<double,3,3> &t3as_2,
	       const Tensor3_antisymmetric<double,3,3> &t3as_3)
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

  test_for_zero(t3as_2(i,j,k)*(t2_2(i,j)*t1_2(k))
		- (t3as_2(0,1,0)*t2_2(0,1)*t1_2(0)
		   + t3as_2(1,1,0)*t2_2(1,1)*t1_2(0)
		   + t3as_2(2,1,0)*t2_2(2,1)*t1_2(0)
		   + t3as_2(0,2,0)*t2_2(0,2)*t1_2(0)
		   + t3as_2(1,2,0)*t2_2(1,2)*t1_2(0)
		   + t3as_2(2,2,0)*t2_2(2,2)*t1_2(0)
		   + t3as_2(0,0,1)*t2_2(0,0)*t1_2(1)
		   + t3as_2(1,0,1)*t2_2(1,0)*t1_2(1)
		   + t3as_2(2,0,1)*t2_2(2,0)*t1_2(1)
		   + t3as_2(0,2,1)*t2_2(0,2)*t1_2(1)
		   + t3as_2(1,2,1)*t2_2(1,2)*t1_2(1)
		   + t3as_2(2,2,1)*t2_2(2,2)*t1_2(1)
		   + t3as_2(0,0,2)*t2_2(0,0)*t1_2(2)
		   + t3as_2(1,0,2)*t2_2(1,0)*t1_2(2)
		   + t3as_2(2,0,2)*t2_2(2,0)*t1_2(2)
		   + t3as_2(0,1,2)*t2_2(0,1)*t1_2(2)
		   + t3as_2(1,1,2)*t2_2(1,1)*t1_2(2)
		   + t3as_2(2,1,2)*t2_2(2,1)*t1_2(2))
		,"T3as(i,j,k)*T3dg(i,j,k)");
  test_for_zero((t2_2(i,j)*t1_2(k)*t3as_2(i,j,k))
		- (t3as_2(0,1,0)*t2_2(0,1)*t1_2(0)
		   + t3as_2(1,1,0)*t2_2(1,1)*t1_2(0)
		   + t3as_2(2,1,0)*t2_2(2,1)*t1_2(0)
		   + t3as_2(0,2,0)*t2_2(0,2)*t1_2(0)
		   + t3as_2(1,2,0)*t2_2(1,2)*t1_2(0)
		   + t3as_2(2,2,0)*t2_2(2,2)*t1_2(0)
		   + t3as_2(0,0,1)*t2_2(0,0)*t1_2(1)
		   + t3as_2(1,0,1)*t2_2(1,0)*t1_2(1)
		   + t3as_2(2,0,1)*t2_2(2,0)*t1_2(1)
		   + t3as_2(0,2,1)*t2_2(0,2)*t1_2(1)
		   + t3as_2(1,2,1)*t2_2(1,2)*t1_2(1)
		   + t3as_2(2,2,1)*t2_2(2,2)*t1_2(1)
		   + t3as_2(0,0,2)*t2_2(0,0)*t1_2(2)
		   + t3as_2(1,0,2)*t2_2(1,0)*t1_2(2)
		   + t3as_2(2,0,2)*t2_2(2,0)*t1_2(2)
		   + t3as_2(0,1,2)*t2_2(0,1)*t1_2(2)
		   + t3as_2(1,1,2)*t2_2(1,1)*t1_2(2)
		   + t3as_2(2,1,2)*t2_2(2,1)*t1_2(2))
		,"T3dg(i,j,k)*T3as(i,j,k)");
}
