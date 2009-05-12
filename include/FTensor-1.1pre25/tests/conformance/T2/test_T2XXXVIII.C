#include <iostream>
#include "../../../FTensor.h"
#include "../test_for_zero.h"
using namespace FTensor;
using namespace std;

void test_T2XXXVIII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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

  /* Checking whether functors work */

  Tensor2<complex<double>,2,2> t1(complex<double>(1,2),complex<double>(3,4),
				  complex<double>(5,6),complex<double>(7,8));

  Index<'a',2> a;
  Index<'b',2> b;

  /* I have this function pointer defined here because I can't seem to
     make the compiler understand what kind of function I have.
     Unless I do casts, which are ugly. It seems to have something to
     do with conj being a templated function. */

  complex<double> (*cj)(const complex<double> &);
  cj=&(conj<double>);

  test_for_zero(conj(t1(0,0))-transform(t1(a,b),static_cast<complex<double>
					(*)(const complex<double> &)>
					(&(conj<double>)))(0,0),
		"transform(T2)(0,0) cast");
  test_for_zero(conj(t1(0,1))-transform(t1(a,b),static_cast<complex<double>
				  (*)(const complex<double> &)>
					(&(conj<double>)))(0,1),
		"transform(T2)(0,1) cast");
  test_for_zero(conj(t1(1,0))-transform(t1(a,b),static_cast<complex<double>
					(*)(const complex<double> &)>
					(&(conj<double>)))(1,0),
		"transform(T2)(1,0) cast");
  test_for_zero(conj(t1(1,1))-transform(t1(a,b),static_cast<complex<double>
					(*)(const complex<double> &)>
					(&(conj<double>)))(1,1),
		"transform(T2)(1,1) cast");
  test_for_zero(conj(t1(0,0))-transform(t1(a,b),cj)(0,0),"transform(T2)(0,0)");
  test_for_zero(conj(t1(0,1))-transform(t1(a,b),cj)(0,1),"transform(T2)(0,1)");
  test_for_zero(conj(t1(1,0))-transform(t1(a,b),cj)(1,0),"transform(T2)(1,0)");
  test_for_zero(conj(t1(1,1))-transform(t1(a,b),cj)(1,1),"transform(T2)(1,1)");

  /* Check plain old conj */

  test_for_zero(conj(t1(0,0))-conj(t1(a,b))(0,0),"conj(T2)(0,0)");
  test_for_zero(conj(t1(0,1))-conj(t1(a,b))(0,1),"conj(T2)(0,1)");
  test_for_zero(conj(t1(1,0))-conj(t1(a,b))(1,0),"conj(T2)(1,0)");
  test_for_zero(conj(t1(1,1))-conj(t1(a,b))(1,1),"conj(T2)(1,1)");

  cout << endl;
}
