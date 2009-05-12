#include <iostream>
#include "../../../FTensor.h"
#include "../test_for_zero.h"
using namespace FTensor;
using namespace std;

extern
void test_T3asI(const int &T, Tensor0<double*> &t0_1,
	       const Tensor0<double*> &t0_2,
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
	       const Tensor3_antisymmetric<double,3,3> &t3as_3);

extern
void test_T3asII(const int &T, Tensor0<double*> &t0_1,
	       const Tensor0<double*> &t0_2,
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
	       const Tensor3_antisymmetric<double,3,3> &t3as_3);

extern
void test_T3asIII(const int &T, Tensor0<double*> &t0_1,
	       const Tensor0<double*> &t0_2,
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
	       const Tensor3_antisymmetric<double,3,3> &t3as_3);

extern
void test_T3asIV(const int &T, Tensor0<double*> &t0_1,
	       const Tensor0<double*> &t0_2,
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
	       const Tensor3_antisymmetric<double,3,3> &t3as_3);

extern
void test_T3asV(const int &T, Tensor0<double*> &t0_1,
	       const Tensor0<double*> &t0_2,
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
	       const Tensor3_antisymmetric<double,3,3> &t3as_3);

extern
void test_T3asVI(const int &T, Tensor0<double*> &t0_1,
	       const Tensor0<double*> &t0_2,
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
	       const Tensor3_antisymmetric<double,3,3> &t3as_3);

extern
void test_T3asVII(const int &T, Tensor0<double*> &t0_1,
	       const Tensor0<double*> &t0_2,
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
	       const Tensor3_antisymmetric<double,3,3> &t3as_3);

extern
void test_T3asVIII(const int &T, Tensor0<double*> &t0_1,
	       const Tensor0<double*> &t0_2,
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
	       const Tensor3_antisymmetric<double,3,3> &t3as_3);

extern
void test_T3asIX(const int &T, Tensor0<double*> &t0_1,
	       const Tensor0<double*> &t0_2,
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
	       const Tensor3_antisymmetric<double,3,3> &t3as_3);

extern
void test_T3asX(const int &T, Tensor0<double*> &t0_1,
	       const Tensor0<double*> &t0_2,
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
	       const Tensor3_antisymmetric<double,3,3> &t3as_3);

extern
void test_T3asXI(const int &T, Tensor0<double*> &t0_1,
	       const Tensor0<double*> &t0_2,
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
	       const Tensor3_antisymmetric<double,3,3> &t3as_3);

extern
void test_T3asXII(const int &T, Tensor0<double*> &t0_1,
	       const Tensor0<double*> &t0_2,
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
	       const Tensor3_antisymmetric<double,3,3> &t3as_3);

extern
void test_T3asXIII(const int &T, Tensor0<double*> &t0_1,
	       const Tensor0<double*> &t0_2,
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
	       const Tensor3_antisymmetric<double,3,3> &t3as_3);

extern
void test_T3asXIV(const int &T, Tensor0<double*> &t0_1,
	       const Tensor0<double*> &t0_2,
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
	       const Tensor3_antisymmetric<double,3,3> &t3as_3);



void test_T3as(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
  test_T3asI(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3as_1,t3as_2,t3as_3);
  test_T3asII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3as_1,t3as_2,t3as_3);
  test_T3asIII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3as_1,t3as_2,t3as_3);
  test_T3asIV(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3as_1,t3as_2,t3as_3);
  test_T3asV(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3as_1,t3as_2,t3as_3);
  test_T3asVI(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3as_1,t3as_2,t3as_3);
  test_T3asVII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3as_1,t3as_2,t3as_3);
  test_T3asVIII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3as_1,t3as_2,t3as_3);
  test_T3asIX(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3as_1,t3as_2,t3as_3);
  test_T3asX(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3as_1,t3as_2,t3as_3);
  test_T3asXI(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3as_1,t3as_2,t3as_3);
  test_T3asXII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3as_1,t3as_2,t3as_3);
  test_T3asXIII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3as_1,t3as_2,t3as_3);
  test_T3asXIV(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3as_1,t3as_2,t3as_3);

  cout << endl;
}
