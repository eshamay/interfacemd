#include <iostream>
#include "../../FTensor.h"
#include "test_for_zero.h"
using namespace FTensor;
using namespace std;

void test_T4R(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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

  Tensor4_ddg<double,3,3> t4ddg_2, t4ddg_3;
  Tensor4_Riemann<double,3> t4R_1, t4R_2, t4R_3;

  cout << "T4_Riemann" << endl;


  t4ddg_2(i,j,l,m)=t2s_2(i,j)*t2s_3(l,m);
  t4ddg_3(i,j,l,m)=t3dg_3(i,j,k)*t3dg_2(l,m,k);

  t4R_1(i,j,k,l)=(t4ddg_2(i,j,k,l)&&t4ddg_2(i,l,k,j));

  test_for_zero(t4R_1(0,0,0,0) - (t4ddg_2(0,0,0,0)-t4ddg_2(0,0,0,0))
      ); test_for_zero(t4R_1(0,0,0,1) - (t4ddg_2(0,0,0,1)-t4ddg_2(0,1,0,0))
      ); test_for_zero(t4R_1(0,0,0,2) - (t4ddg_2(0,0,0,2)-t4ddg_2(0,2,0,0))
      ); test_for_zero(t4R_1(0,0,1,0) - (t4ddg_2(0,0,1,0)-t4ddg_2(0,0,1,0))
      ); test_for_zero(t4R_1(0,0,1,1) - (t4ddg_2(0,0,1,1)-t4ddg_2(0,1,1,0))
      ); test_for_zero(t4R_1(0,0,1,2) - (t4ddg_2(0,0,1,2)-t4ddg_2(0,2,1,0))
      ); test_for_zero(t4R_1(0,0,2,0) - (t4ddg_2(0,0,2,0)-t4ddg_2(0,0,2,0))
      ); test_for_zero(t4R_1(0,0,2,1) - (t4ddg_2(0,0,2,1)-t4ddg_2(0,1,2,0))
      ); test_for_zero(t4R_1(0,0,2,2) - (t4ddg_2(0,0,2,2)-t4ddg_2(0,2,2,0))
      ); test_for_zero(t4R_1(0,1,0,0) - (t4ddg_2(0,1,0,0)-t4ddg_2(0,0,0,1))
      ); test_for_zero(t4R_1(0,1,0,1) - (t4ddg_2(0,1,0,1)-t4ddg_2(0,1,0,1))
      ); test_for_zero(t4R_1(0,1,0,2) - (t4ddg_2(0,1,0,2)-t4ddg_2(0,2,0,1))
      ); test_for_zero(t4R_1(0,1,1,0) - (t4ddg_2(0,1,1,0)-t4ddg_2(0,0,1,1))
      ); test_for_zero(t4R_1(0,1,1,1) - (t4ddg_2(0,1,1,1)-t4ddg_2(0,1,1,1))
      ); test_for_zero(t4R_1(0,1,1,2) - (t4ddg_2(0,1,1,2)-t4ddg_2(0,2,1,1))
      ); test_for_zero(t4R_1(0,1,2,0) - (t4ddg_2(0,1,2,0)-t4ddg_2(0,0,2,1))
      ); test_for_zero(t4R_1(0,1,2,1) - (t4ddg_2(0,1,2,1)-t4ddg_2(0,1,2,1))
      ); test_for_zero(t4R_1(0,1,2,2) - (t4ddg_2(0,1,2,2)-t4ddg_2(0,2,2,1))
      ); test_for_zero(t4R_1(0,2,0,0) - (t4ddg_2(0,2,0,0)-t4ddg_2(0,0,0,2))
      ); test_for_zero(t4R_1(0,2,0,1) - (t4ddg_2(0,2,0,1)-t4ddg_2(0,1,0,2))
      ); test_for_zero(t4R_1(0,2,0,2) - (t4ddg_2(0,2,0,2)-t4ddg_2(0,2,0,2))
      ); test_for_zero(t4R_1(0,2,1,0) - (t4ddg_2(0,2,1,0)-t4ddg_2(0,0,1,2))
      ); test_for_zero(t4R_1(0,2,1,1) - (t4ddg_2(0,2,1,1)-t4ddg_2(0,1,1,2))
      ); test_for_zero(t4R_1(0,2,1,2) - (t4ddg_2(0,2,1,2)-t4ddg_2(0,2,1,2))
      ); test_for_zero(t4R_1(0,2,2,0) - (t4ddg_2(0,2,2,0)-t4ddg_2(0,0,2,2))
      ); test_for_zero(t4R_1(0,2,2,1) - (t4ddg_2(0,2,2,1)-t4ddg_2(0,1,2,2))
      ); test_for_zero(t4R_1(0,2,2,2) - (t4ddg_2(0,2,2,2)-t4ddg_2(0,2,2,2))
      ); test_for_zero(t4R_1(1,0,0,0) - (t4ddg_2(1,0,0,0)-t4ddg_2(1,0,0,0))
      ); test_for_zero(t4R_1(1,0,0,1) - (t4ddg_2(1,0,0,1)-t4ddg_2(1,1,0,0))
      ); test_for_zero(t4R_1(1,0,0,2) - (t4ddg_2(1,0,0,2)-t4ddg_2(1,2,0,0))
      ); test_for_zero(t4R_1(1,0,1,0) - (t4ddg_2(1,0,1,0)-t4ddg_2(1,0,1,0))
      ); test_for_zero(t4R_1(1,0,1,1) - (t4ddg_2(1,0,1,1)-t4ddg_2(1,1,1,0))
      ); test_for_zero(t4R_1(1,0,1,2) - (t4ddg_2(1,0,1,2)-t4ddg_2(1,2,1,0))
      ); test_for_zero(t4R_1(1,0,2,0) - (t4ddg_2(1,0,2,0)-t4ddg_2(1,0,2,0))
      ); test_for_zero(t4R_1(1,0,2,1) - (t4ddg_2(1,0,2,1)-t4ddg_2(1,1,2,0))
      ); test_for_zero(t4R_1(1,0,2,2) - (t4ddg_2(1,0,2,2)-t4ddg_2(1,2,2,0))
      ); test_for_zero(t4R_1(1,1,0,0) - (t4ddg_2(1,1,0,0)-t4ddg_2(1,0,0,1))
      ); test_for_zero(t4R_1(1,1,0,1) - (t4ddg_2(1,1,0,1)-t4ddg_2(1,1,0,1))
      ); test_for_zero(t4R_1(1,1,0,2) - (t4ddg_2(1,1,0,2)-t4ddg_2(1,2,0,1))
      ); test_for_zero(t4R_1(1,1,1,0) - (t4ddg_2(1,1,1,0)-t4ddg_2(1,0,1,1))
      ); test_for_zero(t4R_1(1,1,1,1) - (t4ddg_2(1,1,1,1)-t4ddg_2(1,1,1,1))
      ); test_for_zero(t4R_1(1,1,1,2) - (t4ddg_2(1,1,1,2)-t4ddg_2(1,2,1,1))
      ); test_for_zero(t4R_1(1,1,2,0) - (t4ddg_2(1,1,2,0)-t4ddg_2(1,0,2,1))
      ); test_for_zero(t4R_1(1,1,2,1) - (t4ddg_2(1,1,2,1)-t4ddg_2(1,1,2,1))
      ); test_for_zero(t4R_1(1,1,2,2) - (t4ddg_2(1,1,2,2)-t4ddg_2(1,2,2,1))
      ); test_for_zero(t4R_1(1,2,0,0) - (t4ddg_2(1,2,0,0)-t4ddg_2(1,0,0,2))
      ); test_for_zero(t4R_1(1,2,0,1) - (t4ddg_2(1,2,0,1)-t4ddg_2(1,1,0,2))
      ); test_for_zero(t4R_1(1,2,0,2) - (t4ddg_2(1,2,0,2)-t4ddg_2(1,2,0,2))
      ); test_for_zero(t4R_1(1,2,1,0) - (t4ddg_2(1,2,1,0)-t4ddg_2(1,0,1,2))
      ); test_for_zero(t4R_1(1,2,1,1) - (t4ddg_2(1,2,1,1)-t4ddg_2(1,1,1,2))
      ); test_for_zero(t4R_1(1,2,1,2) - (t4ddg_2(1,2,1,2)-t4ddg_2(1,2,1,2))
      ); test_for_zero(t4R_1(1,2,2,0) - (t4ddg_2(1,2,2,0)-t4ddg_2(1,0,2,2))
      ); test_for_zero(t4R_1(1,2,2,1) - (t4ddg_2(1,2,2,1)-t4ddg_2(1,1,2,2))
      ); test_for_zero(t4R_1(1,2,2,2) - (t4ddg_2(1,2,2,2)-t4ddg_2(1,2,2,2))
      ); test_for_zero(t4R_1(2,0,0,0) - (t4ddg_2(2,0,0,0)-t4ddg_2(2,0,0,0))
      ); test_for_zero(t4R_1(2,0,0,1) - (t4ddg_2(2,0,0,1)-t4ddg_2(2,1,0,0))
      ); test_for_zero(t4R_1(2,0,0,2) - (t4ddg_2(2,0,0,2)-t4ddg_2(2,2,0,0))
      ); test_for_zero(t4R_1(2,0,1,0) - (t4ddg_2(2,0,1,0)-t4ddg_2(2,0,1,0))
      ); test_for_zero(t4R_1(2,0,1,1) - (t4ddg_2(2,0,1,1)-t4ddg_2(2,1,1,0))
      ); test_for_zero(t4R_1(2,0,1,2) - (t4ddg_2(2,0,1,2)-t4ddg_2(2,2,1,0))
      ); test_for_zero(t4R_1(2,0,2,0) - (t4ddg_2(2,0,2,0)-t4ddg_2(2,0,2,0))
      ); test_for_zero(t4R_1(2,0,2,1) - (t4ddg_2(2,0,2,1)-t4ddg_2(2,1,2,0))
      ); test_for_zero(t4R_1(2,0,2,2) - (t4ddg_2(2,0,2,2)-t4ddg_2(2,2,2,0))
      ); test_for_zero(t4R_1(2,1,0,0) - (t4ddg_2(2,1,0,0)-t4ddg_2(2,0,0,1))
      ); test_for_zero(t4R_1(2,1,0,1) - (t4ddg_2(2,1,0,1)-t4ddg_2(2,1,0,1))
      ); test_for_zero(t4R_1(2,1,0,2) - (t4ddg_2(2,1,0,2)-t4ddg_2(2,2,0,1))
      ); test_for_zero(t4R_1(2,1,1,0) - (t4ddg_2(2,1,1,0)-t4ddg_2(2,0,1,1))
      ); test_for_zero(t4R_1(2,1,1,1) - (t4ddg_2(2,1,1,1)-t4ddg_2(2,1,1,1))
      ); test_for_zero(t4R_1(2,1,1,2) - (t4ddg_2(2,1,1,2)-t4ddg_2(2,2,1,1))
      ); test_for_zero(t4R_1(2,1,2,0) - (t4ddg_2(2,1,2,0)-t4ddg_2(2,0,2,1))
      ); test_for_zero(t4R_1(2,1,2,1) - (t4ddg_2(2,1,2,1)-t4ddg_2(2,1,2,1))
      ); test_for_zero(t4R_1(2,1,2,2) - (t4ddg_2(2,1,2,2)-t4ddg_2(2,2,2,1))
      ); test_for_zero(t4R_1(2,2,0,0) - (t4ddg_2(2,2,0,0)-t4ddg_2(2,0,0,2))
      ); test_for_zero(t4R_1(2,2,0,1) - (t4ddg_2(2,2,0,1)-t4ddg_2(2,1,0,2))
      ); test_for_zero(t4R_1(2,2,0,2) - (t4ddg_2(2,2,0,2)-t4ddg_2(2,2,0,2))
      ); test_for_zero(t4R_1(2,2,1,0) - (t4ddg_2(2,2,1,0)-t4ddg_2(2,0,1,2))
      ); test_for_zero(t4R_1(2,2,1,1) - (t4ddg_2(2,2,1,1)-t4ddg_2(2,1,1,2))
      ); test_for_zero(t4R_1(2,2,1,2) - (t4ddg_2(2,2,1,2)-t4ddg_2(2,2,1,2))
      ); test_for_zero(t4R_1(2,2,2,0) - (t4ddg_2(2,2,2,0)-t4ddg_2(2,0,2,2))
      ); test_for_zero(t4R_1(2,2,2,1) - (t4ddg_2(2,2,2,1)-t4ddg_2(2,1,2,2))
      ); test_for_zero(t4R_1(2,2,2,2) - (t4ddg_2(2,2,2,2)-t4ddg_2(2,2,2,2)));

  t4R_2(i,j,k,l)=(t4ddg_3(i,j,k,l)&&t4ddg_3(i,l,k,j));

  t4R_3(i,j,k,l)=t4R_1(i,j,k,l)+t4R_2(i,j,k,l);
  test_for_zero(t4R_3(0,0,0,0) - (t4R_1(0,0,0,0)+t4R_2(0,0,0,0))
      ); test_for_zero(t4R_3(0,0,0,1) - (t4R_1(0,0,0,1)+t4R_2(0,0,0,1))
      ); test_for_zero(t4R_3(0,0,0,2) - (t4R_1(0,0,0,2)+t4R_2(0,0,0,2))
      ); test_for_zero(t4R_3(0,0,1,0) - (t4R_1(0,0,1,0)+t4R_2(0,0,1,0))
      ); test_for_zero(t4R_3(0,0,1,1) - (t4R_1(0,0,1,1)+t4R_2(0,0,1,1))
      ); test_for_zero(t4R_3(0,0,1,2) - (t4R_1(0,0,1,2)+t4R_2(0,0,1,2))
      ); test_for_zero(t4R_3(0,0,2,0) - (t4R_1(0,0,2,0)+t4R_2(0,0,2,0))
      ); test_for_zero(t4R_3(0,0,2,1) - (t4R_1(0,0,2,1)+t4R_2(0,0,2,1))
      ); test_for_zero(t4R_3(0,0,2,2) - (t4R_1(0,0,2,2)+t4R_2(0,0,2,2))
      ); test_for_zero(t4R_3(0,1,0,0) - (t4R_1(0,1,0,0)+t4R_2(0,1,0,0))
      ); test_for_zero(t4R_3(0,1,0,1) - (t4R_1(0,1,0,1)+t4R_2(0,1,0,1))
      ); test_for_zero(t4R_3(0,1,0,2) - (t4R_1(0,1,0,2)+t4R_2(0,1,0,2))
      ); test_for_zero(t4R_3(0,1,1,0) - (t4R_1(0,1,1,0)+t4R_2(0,1,1,0))
      ); test_for_zero(t4R_3(0,1,1,1) - (t4R_1(0,1,1,1)+t4R_2(0,1,1,1))
      ); test_for_zero(t4R_3(0,1,1,2) - (t4R_1(0,1,1,2)+t4R_2(0,1,1,2))
      ); test_for_zero(t4R_3(0,1,2,0) - (t4R_1(0,1,2,0)+t4R_2(0,1,2,0))
      ); test_for_zero(t4R_3(0,1,2,1) - (t4R_1(0,1,2,1)+t4R_2(0,1,2,1))
      ); test_for_zero(t4R_3(0,1,2,2) - (t4R_1(0,1,2,2)+t4R_2(0,1,2,2))
      ); test_for_zero(t4R_3(0,2,0,0) - (t4R_1(0,2,0,0)+t4R_2(0,2,0,0))
      ); test_for_zero(t4R_3(0,2,0,1) - (t4R_1(0,2,0,1)+t4R_2(0,2,0,1))
      ); test_for_zero(t4R_3(0,2,0,2) - (t4R_1(0,2,0,2)+t4R_2(0,2,0,2))
      ); test_for_zero(t4R_3(0,2,1,0) - (t4R_1(0,2,1,0)+t4R_2(0,2,1,0))
      ); test_for_zero(t4R_3(0,2,1,1) - (t4R_1(0,2,1,1)+t4R_2(0,2,1,1))
      ); test_for_zero(t4R_3(0,2,1,2) - (t4R_1(0,2,1,2)+t4R_2(0,2,1,2))
      ); test_for_zero(t4R_3(0,2,2,0) - (t4R_1(0,2,2,0)+t4R_2(0,2,2,0))
      ); test_for_zero(t4R_3(0,2,2,1) - (t4R_1(0,2,2,1)+t4R_2(0,2,2,1))
      ); test_for_zero(t4R_3(0,2,2,2) - (t4R_1(0,2,2,2)+t4R_2(0,2,2,2))
      ); test_for_zero(t4R_3(1,0,0,0) - (t4R_1(1,0,0,0)+t4R_2(1,0,0,0))
      ); test_for_zero(t4R_3(1,0,0,1) - (t4R_1(1,0,0,1)+t4R_2(1,0,0,1))
      ); test_for_zero(t4R_3(1,0,0,2) - (t4R_1(1,0,0,2)+t4R_2(1,0,0,2))
      ); test_for_zero(t4R_3(1,0,1,0) - (t4R_1(1,0,1,0)+t4R_2(1,0,1,0))
      ); test_for_zero(t4R_3(1,0,1,1) - (t4R_1(1,0,1,1)+t4R_2(1,0,1,1))
      ); test_for_zero(t4R_3(1,0,1,2) - (t4R_1(1,0,1,2)+t4R_2(1,0,1,2))
      ); test_for_zero(t4R_3(1,0,2,0) - (t4R_1(1,0,2,0)+t4R_2(1,0,2,0))
      ); test_for_zero(t4R_3(1,0,2,1) - (t4R_1(1,0,2,1)+t4R_2(1,0,2,1))
      ); test_for_zero(t4R_3(1,0,2,2) - (t4R_1(1,0,2,2)+t4R_2(1,0,2,2))
      ); test_for_zero(t4R_3(1,1,0,0) - (t4R_1(1,1,0,0)+t4R_2(1,1,0,0))
      ); test_for_zero(t4R_3(1,1,0,1) - (t4R_1(1,1,0,1)+t4R_2(1,1,0,1))
      ); test_for_zero(t4R_3(1,1,0,2) - (t4R_1(1,1,0,2)+t4R_2(1,1,0,2))
      ); test_for_zero(t4R_3(1,1,1,0) - (t4R_1(1,1,1,0)+t4R_2(1,1,1,0))
      ); test_for_zero(t4R_3(1,1,1,1) - (t4R_1(1,1,1,1)+t4R_2(1,1,1,1))
      ); test_for_zero(t4R_3(1,1,1,2) - (t4R_1(1,1,1,2)+t4R_2(1,1,1,2))
      ); test_for_zero(t4R_3(1,1,2,0) - (t4R_1(1,1,2,0)+t4R_2(1,1,2,0))
      ); test_for_zero(t4R_3(1,1,2,1) - (t4R_1(1,1,2,1)+t4R_2(1,1,2,1))
      ); test_for_zero(t4R_3(1,1,2,2) - (t4R_1(1,1,2,2)+t4R_2(1,1,2,2))
      ); test_for_zero(t4R_3(1,2,0,0) - (t4R_1(1,2,0,0)+t4R_2(1,2,0,0))
      ); test_for_zero(t4R_3(1,2,0,1) - (t4R_1(1,2,0,1)+t4R_2(1,2,0,1))
      ); test_for_zero(t4R_3(1,2,0,2) - (t4R_1(1,2,0,2)+t4R_2(1,2,0,2))
      ); test_for_zero(t4R_3(1,2,1,0) - (t4R_1(1,2,1,0)+t4R_2(1,2,1,0))
      ); test_for_zero(t4R_3(1,2,1,1) - (t4R_1(1,2,1,1)+t4R_2(1,2,1,1))
      ); test_for_zero(t4R_3(1,2,1,2) - (t4R_1(1,2,1,2)+t4R_2(1,2,1,2))
      ); test_for_zero(t4R_3(1,2,2,0) - (t4R_1(1,2,2,0)+t4R_2(1,2,2,0))
      ); test_for_zero(t4R_3(1,2,2,1) - (t4R_1(1,2,2,1)+t4R_2(1,2,2,1))
      ); test_for_zero(t4R_3(1,2,2,2) - (t4R_1(1,2,2,2)+t4R_2(1,2,2,2))
      ); test_for_zero(t4R_3(2,0,0,0) - (t4R_1(2,0,0,0)+t4R_2(2,0,0,0))
      ); test_for_zero(t4R_3(2,0,0,1) - (t4R_1(2,0,0,1)+t4R_2(2,0,0,1))
      ); test_for_zero(t4R_3(2,0,0,2) - (t4R_1(2,0,0,2)+t4R_2(2,0,0,2))
      ); test_for_zero(t4R_3(2,0,1,0) - (t4R_1(2,0,1,0)+t4R_2(2,0,1,0))
      ); test_for_zero(t4R_3(2,0,1,1) - (t4R_1(2,0,1,1)+t4R_2(2,0,1,1))
      ); test_for_zero(t4R_3(2,0,1,2) - (t4R_1(2,0,1,2)+t4R_2(2,0,1,2))
      ); test_for_zero(t4R_3(2,0,2,0) - (t4R_1(2,0,2,0)+t4R_2(2,0,2,0))
      ); test_for_zero(t4R_3(2,0,2,1) - (t4R_1(2,0,2,1)+t4R_2(2,0,2,1))
      ); test_for_zero(t4R_3(2,0,2,2) - (t4R_1(2,0,2,2)+t4R_2(2,0,2,2))
      ); test_for_zero(t4R_3(2,1,0,0) - (t4R_1(2,1,0,0)+t4R_2(2,1,0,0))
      ); test_for_zero(t4R_3(2,1,0,1) - (t4R_1(2,1,0,1)+t4R_2(2,1,0,1))
      ); test_for_zero(t4R_3(2,1,0,2) - (t4R_1(2,1,0,2)+t4R_2(2,1,0,2))
      ); test_for_zero(t4R_3(2,1,1,0) - (t4R_1(2,1,1,0)+t4R_2(2,1,1,0))
      ); test_for_zero(t4R_3(2,1,1,1) - (t4R_1(2,1,1,1)+t4R_2(2,1,1,1))
      ); test_for_zero(t4R_3(2,1,1,2) - (t4R_1(2,1,1,2)+t4R_2(2,1,1,2))
      ); test_for_zero(t4R_3(2,1,2,0) - (t4R_1(2,1,2,0)+t4R_2(2,1,2,0))
      ); test_for_zero(t4R_3(2,1,2,1) - (t4R_1(2,1,2,1)+t4R_2(2,1,2,1))
      ); test_for_zero(t4R_3(2,1,2,2) - (t4R_1(2,1,2,2)+t4R_2(2,1,2,2))
      ); test_for_zero(t4R_3(2,2,0,0) - (t4R_1(2,2,0,0)+t4R_2(2,2,0,0))
      ); test_for_zero(t4R_3(2,2,0,1) - (t4R_1(2,2,0,1)+t4R_2(2,2,0,1))
      ); test_for_zero(t4R_3(2,2,0,2) - (t4R_1(2,2,0,2)+t4R_2(2,2,0,2))
      ); test_for_zero(t4R_3(2,2,1,0) - (t4R_1(2,2,1,0)+t4R_2(2,2,1,0))
      ); test_for_zero(t4R_3(2,2,1,1) - (t4R_1(2,2,1,1)+t4R_2(2,2,1,1))
      ); test_for_zero(t4R_3(2,2,1,2) - (t4R_1(2,2,1,2)+t4R_2(2,2,1,2))
      ); test_for_zero(t4R_3(2,2,2,0) - (t4R_1(2,2,2,0)+t4R_2(2,2,2,0))
      ); test_for_zero(t4R_3(2,2,2,1) - (t4R_1(2,2,2,1)+t4R_2(2,2,2,1))
      ); test_for_zero(t4R_3(2,2,2,2) - (t4R_1(2,2,2,2)+t4R_2(2,2,2,2)));


  t4R_3(i,j,k,l)=t4R_1(i,j,k,l)-t4R_2(i,j,k,l);
  test_for_zero(t4R_3(0,0,0,0) - (t4R_1(0,0,0,0)-t4R_2(0,0,0,0))
      ); test_for_zero(t4R_3(0,0,0,1) - (t4R_1(0,0,0,1)-t4R_2(0,0,0,1))
      ); test_for_zero(t4R_3(0,0,0,2) - (t4R_1(0,0,0,2)-t4R_2(0,0,0,2))
      ); test_for_zero(t4R_3(0,0,1,0) - (t4R_1(0,0,1,0)-t4R_2(0,0,1,0))
      ); test_for_zero(t4R_3(0,0,1,1) - (t4R_1(0,0,1,1)-t4R_2(0,0,1,1))
      ); test_for_zero(t4R_3(0,0,1,2) - (t4R_1(0,0,1,2)-t4R_2(0,0,1,2))
      ); test_for_zero(t4R_3(0,0,2,0) - (t4R_1(0,0,2,0)-t4R_2(0,0,2,0))
      ); test_for_zero(t4R_3(0,0,2,1) - (t4R_1(0,0,2,1)-t4R_2(0,0,2,1))
      ); test_for_zero(t4R_3(0,0,2,2) - (t4R_1(0,0,2,2)-t4R_2(0,0,2,2))
      ); test_for_zero(t4R_3(0,1,0,0) - (t4R_1(0,1,0,0)-t4R_2(0,1,0,0))
      ); test_for_zero(t4R_3(0,1,0,1) - (t4R_1(0,1,0,1)-t4R_2(0,1,0,1))
      ); test_for_zero(t4R_3(0,1,0,2) - (t4R_1(0,1,0,2)-t4R_2(0,1,0,2))
      ); test_for_zero(t4R_3(0,1,1,0) - (t4R_1(0,1,1,0)-t4R_2(0,1,1,0))
      ); test_for_zero(t4R_3(0,1,1,1) - (t4R_1(0,1,1,1)-t4R_2(0,1,1,1))
      ); test_for_zero(t4R_3(0,1,1,2) - (t4R_1(0,1,1,2)-t4R_2(0,1,1,2))
      ); test_for_zero(t4R_3(0,1,2,0) - (t4R_1(0,1,2,0)-t4R_2(0,1,2,0))
      ); test_for_zero(t4R_3(0,1,2,1) - (t4R_1(0,1,2,1)-t4R_2(0,1,2,1))
      ); test_for_zero(t4R_3(0,1,2,2) - (t4R_1(0,1,2,2)-t4R_2(0,1,2,2))
      ); test_for_zero(t4R_3(0,2,0,0) - (t4R_1(0,2,0,0)-t4R_2(0,2,0,0))
      ); test_for_zero(t4R_3(0,2,0,1) - (t4R_1(0,2,0,1)-t4R_2(0,2,0,1))
      ); test_for_zero(t4R_3(0,2,0,2) - (t4R_1(0,2,0,2)-t4R_2(0,2,0,2))
      ); test_for_zero(t4R_3(0,2,1,0) - (t4R_1(0,2,1,0)-t4R_2(0,2,1,0))
      ); test_for_zero(t4R_3(0,2,1,1) - (t4R_1(0,2,1,1)-t4R_2(0,2,1,1))
      ); test_for_zero(t4R_3(0,2,1,2) - (t4R_1(0,2,1,2)-t4R_2(0,2,1,2))
      ); test_for_zero(t4R_3(0,2,2,0) - (t4R_1(0,2,2,0)-t4R_2(0,2,2,0))
      ); test_for_zero(t4R_3(0,2,2,1) - (t4R_1(0,2,2,1)-t4R_2(0,2,2,1))
      ); test_for_zero(t4R_3(0,2,2,2) - (t4R_1(0,2,2,2)-t4R_2(0,2,2,2))
      ); test_for_zero(t4R_3(1,0,0,0) - (t4R_1(1,0,0,0)-t4R_2(1,0,0,0))
      ); test_for_zero(t4R_3(1,0,0,1) - (t4R_1(1,0,0,1)-t4R_2(1,0,0,1))
      ); test_for_zero(t4R_3(1,0,0,2) - (t4R_1(1,0,0,2)-t4R_2(1,0,0,2))
      ); test_for_zero(t4R_3(1,0,1,0) - (t4R_1(1,0,1,0)-t4R_2(1,0,1,0))
      ); test_for_zero(t4R_3(1,0,1,1) - (t4R_1(1,0,1,1)-t4R_2(1,0,1,1))
      ); test_for_zero(t4R_3(1,0,1,2) - (t4R_1(1,0,1,2)-t4R_2(1,0,1,2))
      ); test_for_zero(t4R_3(1,0,2,0) - (t4R_1(1,0,2,0)-t4R_2(1,0,2,0))
      ); test_for_zero(t4R_3(1,0,2,1) - (t4R_1(1,0,2,1)-t4R_2(1,0,2,1))
      ); test_for_zero(t4R_3(1,0,2,2) - (t4R_1(1,0,2,2)-t4R_2(1,0,2,2))
      ); test_for_zero(t4R_3(1,1,0,0) - (t4R_1(1,1,0,0)-t4R_2(1,1,0,0))
      ); test_for_zero(t4R_3(1,1,0,1) - (t4R_1(1,1,0,1)-t4R_2(1,1,0,1))
      ); test_for_zero(t4R_3(1,1,0,2) - (t4R_1(1,1,0,2)-t4R_2(1,1,0,2))
      ); test_for_zero(t4R_3(1,1,1,0) - (t4R_1(1,1,1,0)-t4R_2(1,1,1,0))
      ); test_for_zero(t4R_3(1,1,1,1) - (t4R_1(1,1,1,1)-t4R_2(1,1,1,1))
      ); test_for_zero(t4R_3(1,1,1,2) - (t4R_1(1,1,1,2)-t4R_2(1,1,1,2))
      ); test_for_zero(t4R_3(1,1,2,0) - (t4R_1(1,1,2,0)-t4R_2(1,1,2,0))
      ); test_for_zero(t4R_3(1,1,2,1) - (t4R_1(1,1,2,1)-t4R_2(1,1,2,1))
      ); test_for_zero(t4R_3(1,1,2,2) - (t4R_1(1,1,2,2)-t4R_2(1,1,2,2))
      ); test_for_zero(t4R_3(1,2,0,0) - (t4R_1(1,2,0,0)-t4R_2(1,2,0,0))
      ); test_for_zero(t4R_3(1,2,0,1) - (t4R_1(1,2,0,1)-t4R_2(1,2,0,1))
      ); test_for_zero(t4R_3(1,2,0,2) - (t4R_1(1,2,0,2)-t4R_2(1,2,0,2))
      ); test_for_zero(t4R_3(1,2,1,0) - (t4R_1(1,2,1,0)-t4R_2(1,2,1,0))
      ); test_for_zero(t4R_3(1,2,1,1) - (t4R_1(1,2,1,1)-t4R_2(1,2,1,1))
      ); test_for_zero(t4R_3(1,2,1,2) - (t4R_1(1,2,1,2)-t4R_2(1,2,1,2))
      ); test_for_zero(t4R_3(1,2,2,0) - (t4R_1(1,2,2,0)-t4R_2(1,2,2,0))
      ); test_for_zero(t4R_3(1,2,2,1) - (t4R_1(1,2,2,1)-t4R_2(1,2,2,1))
      ); test_for_zero(t4R_3(1,2,2,2) - (t4R_1(1,2,2,2)-t4R_2(1,2,2,2))
      ); test_for_zero(t4R_3(2,0,0,0) - (t4R_1(2,0,0,0)-t4R_2(2,0,0,0))
      ); test_for_zero(t4R_3(2,0,0,1) - (t4R_1(2,0,0,1)-t4R_2(2,0,0,1))
      ); test_for_zero(t4R_3(2,0,0,2) - (t4R_1(2,0,0,2)-t4R_2(2,0,0,2))
      ); test_for_zero(t4R_3(2,0,1,0) - (t4R_1(2,0,1,0)-t4R_2(2,0,1,0))
      ); test_for_zero(t4R_3(2,0,1,1) - (t4R_1(2,0,1,1)-t4R_2(2,0,1,1))
      ); test_for_zero(t4R_3(2,0,1,2) - (t4R_1(2,0,1,2)-t4R_2(2,0,1,2))
      ); test_for_zero(t4R_3(2,0,2,0) - (t4R_1(2,0,2,0)-t4R_2(2,0,2,0))
      ); test_for_zero(t4R_3(2,0,2,1) - (t4R_1(2,0,2,1)-t4R_2(2,0,2,1))
      ); test_for_zero(t4R_3(2,0,2,2) - (t4R_1(2,0,2,2)-t4R_2(2,0,2,2))
      ); test_for_zero(t4R_3(2,1,0,0) - (t4R_1(2,1,0,0)-t4R_2(2,1,0,0))
      ); test_for_zero(t4R_3(2,1,0,1) - (t4R_1(2,1,0,1)-t4R_2(2,1,0,1))
      ); test_for_zero(t4R_3(2,1,0,2) - (t4R_1(2,1,0,2)-t4R_2(2,1,0,2))
      ); test_for_zero(t4R_3(2,1,1,0) - (t4R_1(2,1,1,0)-t4R_2(2,1,1,0))
      ); test_for_zero(t4R_3(2,1,1,1) - (t4R_1(2,1,1,1)-t4R_2(2,1,1,1))
      ); test_for_zero(t4R_3(2,1,1,2) - (t4R_1(2,1,1,2)-t4R_2(2,1,1,2))
      ); test_for_zero(t4R_3(2,1,2,0) - (t4R_1(2,1,2,0)-t4R_2(2,1,2,0))
      ); test_for_zero(t4R_3(2,1,2,1) - (t4R_1(2,1,2,1)-t4R_2(2,1,2,1))
      ); test_for_zero(t4R_3(2,1,2,2) - (t4R_1(2,1,2,2)-t4R_2(2,1,2,2))
      ); test_for_zero(t4R_3(2,2,0,0) - (t4R_1(2,2,0,0)-t4R_2(2,2,0,0))
      ); test_for_zero(t4R_3(2,2,0,1) - (t4R_1(2,2,0,1)-t4R_2(2,2,0,1))
      ); test_for_zero(t4R_3(2,2,0,2) - (t4R_1(2,2,0,2)-t4R_2(2,2,0,2))
      ); test_for_zero(t4R_3(2,2,1,0) - (t4R_1(2,2,1,0)-t4R_2(2,2,1,0))
      ); test_for_zero(t4R_3(2,2,1,1) - (t4R_1(2,2,1,1)-t4R_2(2,2,1,1))
      ); test_for_zero(t4R_3(2,2,1,2) - (t4R_1(2,2,1,2)-t4R_2(2,2,1,2))
      ); test_for_zero(t4R_3(2,2,2,0) - (t4R_1(2,2,2,0)-t4R_2(2,2,2,0))
      ); test_for_zero(t4R_3(2,2,2,1) - (t4R_1(2,2,2,1)-t4R_2(2,2,2,1))
      ); test_for_zero(t4R_3(2,2,2,2) - (t4R_1(2,2,2,2)-t4R_2(2,2,2,2)));

  t3as_1(i,j,k)=t4R_1(i,j,k,l)*t1_2(l);
  t3as_1(i,j,k)=t1_2(l)*t4R_1(i,j,k,l);
  t3as_1(i,j,k)=t4R_1(i,j,l,k)*t1_2(l);
  t3as_1(i,j,k)=t1_2(l)*t4R_1(i,j,l,k);
  t3as_1(i,j,k)=t4R_1(i,l,j,k)*t1_2(l);
  t3as_1(i,j,k)=t1_2(l)*t4R_1(i,l,j,k);
  t3as_1(i,j,k)=t4R_1(l,i,j,k)*t1_2(l);
  t3as_1(i,j,k)=t1_2(l)*t4R_1(l,i,j,k);

  cout << '\n';

  test_for_zero(t4R_1(i,j,k,l)*t4ddg_2(i,j,k,l)
      ); test_for_zero(t4ddg_2(i,j,k,l)*t4R_1(i,j,k,l)
      ); test_for_zero(t4R_1(i,j,k,l)*t4ddg_2(i,k,j,l)
      ); test_for_zero(t4ddg_2(i,k,j,l)*t4R_1(i,j,k,l)
      ); test_for_zero(t4R_3(i,j,k,l)*t4ddg_2(i,j,k,l)
      ); test_for_zero(t4ddg_2(i,j,k,l)*t4R_3(i,j,k,l)
      ); test_for_zero(t4R_3(i,j,k,l)*t4ddg_2(i,k,j,l)
      ); test_for_zero(t4ddg_2(i,k,j,l)*t4R_3(i,j,k,l));

  cout << '\n';

  t2s_1(j,l)=t4R_1(i,j,k,l)*t2s_2(i,k);
  t2s_1(j,l)=t2s_2(i,k)*t4R_1(i,j,k,l);
}
