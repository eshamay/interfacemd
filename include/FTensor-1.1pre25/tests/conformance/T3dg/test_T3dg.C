#include <iostream>
#include "../../../FTensor.h"
#include "../test_for_zero.h"
using namespace FTensor;
using namespace std;

extern
void test_T3dgI(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgIII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgIV(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgV(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgVI(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgVII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgVIII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgIX(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgXI(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgXII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgXIII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgXIV(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgXV(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgXVI(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgXVII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgXVIII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgXIX(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgXX(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgXXI(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgXXII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgXXIII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgXXIV(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgXXV(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgXXVI(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgXXVII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgXXVIII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgXXIX(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgXXX(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgXXXI(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgXXXII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgXXXIII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgXXXIV(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgXXXV(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgCII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);

extern
void test_T3dgCIII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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
	       const Tensor3_christof<double,3,3> &t3ch_3);


void test_T3dg(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
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

  test_T3dgI(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgIII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgIV(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgV(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgVI(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgVII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgVIII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgIX(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgX(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgXI(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgXII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgXIII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgXIV(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgXV(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgXVI(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgXVII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgXVIII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgXIX(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgXX(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgXXI(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgXXII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgXXIII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgXXIV(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgXXV(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgXXVI(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgXXVII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgXXVIII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgXXIX(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgXXX(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgXXXI(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgXXXII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgXXXIII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,
		  t3dg_1,t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgXXXIV(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgXXXV(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);

  test_T3dgCII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);
  test_T3dgCIII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3,t3dg_1,
	     t3dg_2,t3dg_3,t3ch_1,t3ch_2,t3ch_3);

  cout << endl;
}

