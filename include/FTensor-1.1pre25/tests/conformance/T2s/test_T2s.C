#include <iostream>
#include "../../../FTensor.h"
#include "../test_for_zero.h"
using namespace FTensor;
using namespace std;

extern
void test_T2sI(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sIII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sIV(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sV(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sVI(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sVII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sVIII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sIX(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sX(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXI(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXIII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXIV(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXV(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXVI(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXVII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXVIII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXIX(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXX(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXXI(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXXII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXXIII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXXIV(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXXV(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXXVI(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXXVII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXXVIII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXXIX(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXXX(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXXXI(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXXXII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXXXIII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXXXIV(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXXXV(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXXXVI(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXXXVII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXXXVIII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXXXIX(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXL(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXLI(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXLII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXLIII(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

extern
void test_T2sXLIV(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3);

void test_T2s(const int &T, Tensor0<double*> &t0_1, const Tensor0<double*> &t0_2,
	      Tensor1<double,3> &t1_1, const Tensor1<double,3> &t1_2,
	      Tensor2<double,3,3> &t2_1, const Tensor2<double,3,3> &t2_2,
	      const Tensor2<double,3,3> &t2_3,
	      Tensor2_symmetric<double,3> &t2s_1,
	      const Tensor2_symmetric<double,3> &t2s_2,
	      const Tensor2_symmetric<double,3> &t2s_3)
{
  test_T2sI(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sIII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sIV(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sV(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sVI(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sVII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sVIII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sIX(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sX(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXI(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXIII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXIV(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXV(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXVI(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXVII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXVIII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXIX(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXX(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXXI(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXXII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXXIII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXXIV(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXXV(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXXVI(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXXVII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXXVIII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXXIX(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXXX(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXXXI(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXXXII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXXXIII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXXXIV(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXXXV(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXXXVI(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXXXVII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXXXVIII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXXXIX(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXL(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXLI(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXLII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXLIII(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);
  test_T2sXLIV(T,t0_1,t0_2,t1_1,t1_2,t2_1,t2_2,t2_3,t2s_1,t2s_2,t2s_3);

  cout << endl;
}
