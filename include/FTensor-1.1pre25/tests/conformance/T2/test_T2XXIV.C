#include <iostream>
#include "../../../FTensor.h"
#include "../test_for_zero.h"
using namespace FTensor;
using namespace std;

void test_T2XXIV(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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

  /* Tensor2*Tensor1 */

  test_for_zero((t2_2(i,j)*t1_2(k))(0,0,0) - t2_2(0,0)*t1_2(0),
		"T2(i,j)*T1(k)(0,0,0)");
  test_for_zero((t2_2(i,j)*t1_2(k))(0,0,1) - t2_2(0,0)*t1_2(1),
		"T2(i,j)*T1(k)(0,0,1)");
  test_for_zero((t2_2(i,j)*t1_2(k))(0,0,2) - t2_2(0,0)*t1_2(2),
		"T2(i,j)*T1(k)(0,0,2)");
  test_for_zero((t2_2(i,j)*t1_2(k))(0,1,0) - t2_2(0,1)*t1_2(0),
		"T2(i,j)*T1(k)(0,1,0)");
  test_for_zero((t2_2(i,j)*t1_2(k))(0,1,1) - t2_2(0,1)*t1_2(1),
		"T2(i,j)*T1(k)(0,1,1)");
  test_for_zero((t2_2(i,j)*t1_2(k))(0,1,2) - t2_2(0,1)*t1_2(2),
		"T2(i,j)*T1(k)(0,1,2)");
  test_for_zero((t2_2(i,j)*t1_2(k))(0,2,0) - t2_2(0,2)*t1_2(0),
		"T2(i,j)*T1(k)(0,2,0)");
  test_for_zero((t2_2(i,j)*t1_2(k))(0,2,1) - t2_2(0,2)*t1_2(1),
		"T2(i,j)*T1(k)(0,2,1)");
  test_for_zero((t2_2(i,j)*t1_2(k))(0,2,2) - t2_2(0,2)*t1_2(2),
		"T2(i,j)*T1(k)(0,2,2)");
  test_for_zero((t2_2(i,j)*t1_2(k))(1,0,0) - t2_2(1,0)*t1_2(0),
		"T2(i,j)*T1(k)(1,0,0)");
  test_for_zero((t2_2(i,j)*t1_2(k))(1,0,1) - t2_2(1,0)*t1_2(1),
		"T2(i,j)*T1(k)(1,0,1)");
  test_for_zero((t2_2(i,j)*t1_2(k))(1,0,2) - t2_2(1,0)*t1_2(2),
		"T2(i,j)*T1(k)(1,0,2)");
  test_for_zero((t2_2(i,j)*t1_2(k))(1,1,0) - t2_2(1,1)*t1_2(0),
		"T2(i,j)*T1(k)(1,1,0)");
  test_for_zero((t2_2(i,j)*t1_2(k))(1,1,1) - t2_2(1,1)*t1_2(1),
		"T2(i,j)*T1(k)(1,1,1)");
  test_for_zero((t2_2(i,j)*t1_2(k))(1,1,2) - t2_2(1,1)*t1_2(2),
		"T2(i,j)*T1(k)(1,1,2)");
  test_for_zero((t2_2(i,j)*t1_2(k))(1,2,0) - t2_2(1,2)*t1_2(0),
		"T2(i,j)*T1(k)(1,2,0)");
  test_for_zero((t2_2(i,j)*t1_2(k))(1,2,1) - t2_2(1,2)*t1_2(1),
		"T2(i,j)*T1(k)(1,2,1)");
  test_for_zero((t2_2(i,j)*t1_2(k))(1,2,2) - t2_2(1,2)*t1_2(2),
		"T2(i,j)*T1(k)(1,2,2)");
  test_for_zero((t2_2(i,j)*t1_2(k))(2,0,0) - t2_2(2,0)*t1_2(0),
		"T2(i,j)*T1(k)(2,0,0)");
  test_for_zero((t2_2(i,j)*t1_2(k))(2,0,1) - t2_2(2,0)*t1_2(1),
		"T2(i,j)*T1(k)(2,0,1)");
  test_for_zero((t2_2(i,j)*t1_2(k))(2,0,2) - t2_2(2,0)*t1_2(2),
		"T2(i,j)*T1(k)(2,0,2)");
  test_for_zero((t2_2(i,j)*t1_2(k))(2,1,0) - t2_2(2,1)*t1_2(0),
		"T2(i,j)*T1(k)(2,1,0)");
  test_for_zero((t2_2(i,j)*t1_2(k))(2,1,1) - t2_2(2,1)*t1_2(1),
		"T2(i,j)*T1(k)(2,1,1)");
  test_for_zero((t2_2(i,j)*t1_2(k))(2,1,2) - t2_2(2,1)*t1_2(2),
		"T2(i,j)*T1(k)(2,1,2)");
  test_for_zero((t2_2(i,j)*t1_2(k))(2,2,0) - t2_2(2,2)*t1_2(0),
		"T2(i,j)*T1(k)(2,2,0)");
  test_for_zero((t2_2(i,j)*t1_2(k))(2,2,1) - t2_2(2,2)*t1_2(1),
		"T2(i,j)*T1(k)(2,2,1)");
  test_for_zero((t2_2(i,j)*t1_2(k))(2,2,2) - t2_2(2,2)*t1_2(2),
		"T2(i,j)*T1(k)(2,2,2)");

  test_for_zero((t1_2(k)*t2_2(i,j))(0,0,0) - t2_2(0,0)*t1_2(0),
		"T1(k)*T2(i,j)(0,0,0)");
  test_for_zero((t1_2(k)*t2_2(i,j))(0,0,1) - t2_2(0,0)*t1_2(1),
		"T1(k)*T2(i,j)(0,0,1)");
  test_for_zero((t1_2(k)*t2_2(i,j))(0,0,2) - t2_2(0,0)*t1_2(2),
		"T1(k)*T2(i,j)(0,0,2)");
  test_for_zero((t1_2(k)*t2_2(i,j))(0,1,0) - t2_2(0,1)*t1_2(0),
		"T1(k)*T2(i,j)(0,1,0)");
  test_for_zero((t1_2(k)*t2_2(i,j))(0,1,1) - t2_2(0,1)*t1_2(1),
		"T1(k)*T2(i,j)(0,1,1)");
  test_for_zero((t1_2(k)*t2_2(i,j))(0,1,2) - t2_2(0,1)*t1_2(2),
		"T1(k)*T2(i,j)(0,1,2)");
  test_for_zero((t1_2(k)*t2_2(i,j))(0,2,0) - t2_2(0,2)*t1_2(0),
		"T1(k)*T2(i,j)(0,2,0)");
  test_for_zero((t1_2(k)*t2_2(i,j))(0,2,1) - t2_2(0,2)*t1_2(1),
		"T1(k)*T2(i,j)(0,2,1)");
  test_for_zero((t1_2(k)*t2_2(i,j))(0,2,2) - t2_2(0,2)*t1_2(2),
		"T1(k)*T2(i,j)(0,2,2)");
  test_for_zero((t1_2(k)*t2_2(i,j))(1,0,0) - t2_2(1,0)*t1_2(0),
		"T1(k)*T2(i,j)(1,0,0)");
  test_for_zero((t1_2(k)*t2_2(i,j))(1,0,1) - t2_2(1,0)*t1_2(1),
		"T1(k)*T2(i,j)(1,0,1)");
  test_for_zero((t1_2(k)*t2_2(i,j))(1,0,2) - t2_2(1,0)*t1_2(2),
		"T1(k)*T2(i,j)(1,0,2)");
  test_for_zero((t1_2(k)*t2_2(i,j))(1,1,0) - t2_2(1,1)*t1_2(0),
		"T1(k)*T2(i,j)(1,1,0)");
  test_for_zero((t1_2(k)*t2_2(i,j))(1,1,1) - t2_2(1,1)*t1_2(1),
		"T1(k)*T2(i,j)(1,1,1)");
  test_for_zero((t1_2(k)*t2_2(i,j))(1,1,2) - t2_2(1,1)*t1_2(2),
		"T1(k)*T2(i,j)(1,1,2)");
  test_for_zero((t1_2(k)*t2_2(i,j))(1,2,0) - t2_2(1,2)*t1_2(0),
		"T1(k)*T2(i,j)(1,2,0)");
  test_for_zero((t1_2(k)*t2_2(i,j))(1,2,1) - t2_2(1,2)*t1_2(1),
		"T1(k)*T2(i,j)(1,2,1)");
  test_for_zero((t1_2(k)*t2_2(i,j))(1,2,2) - t2_2(1,2)*t1_2(2),
		"T1(k)*T2(i,j)(1,2,2)");
  test_for_zero((t1_2(k)*t2_2(i,j))(2,0,0) - t2_2(2,0)*t1_2(0),
		"T1(k)*T2(i,j)(2,0,0)");
  test_for_zero((t1_2(k)*t2_2(i,j))(2,0,1) - t2_2(2,0)*t1_2(1),
		"T1(k)*T2(i,j)(2,0,1)");
  test_for_zero((t1_2(k)*t2_2(i,j))(2,0,2) - t2_2(2,0)*t1_2(2),
		"T1(k)*T2(i,j)(2,0,2)");
  test_for_zero((t1_2(k)*t2_2(i,j))(2,1,0) - t2_2(2,1)*t1_2(0),
		"T1(k)*T2(i,j)(2,1,0)");
  test_for_zero((t1_2(k)*t2_2(i,j))(2,1,1) - t2_2(2,1)*t1_2(1),
		"T1(k)*T2(i,j)(2,1,1)");
  test_for_zero((t1_2(k)*t2_2(i,j))(2,1,2) - t2_2(2,1)*t1_2(2),
		"T1(k)*T2(i,j)(2,1,2)");
  test_for_zero((t1_2(k)*t2_2(i,j))(2,2,0) - t2_2(2,2)*t1_2(0),
		"T1(k)*T2(i,j)(2,2,0)");
  test_for_zero((t1_2(k)*t2_2(i,j))(2,2,1) - t2_2(2,2)*t1_2(1),
		"T1(k)*T2(i,j)(2,2,1)");
  test_for_zero((t1_2(k)*t2_2(i,j))(2,2,2) - t2_2(2,2)*t1_2(2),
		"T1(k)*T2(i,j)(2,2,2)");
}
