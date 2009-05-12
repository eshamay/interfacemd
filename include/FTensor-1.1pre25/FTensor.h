/* Include file for the Fast Tensor classes (FTensor).  Everything is
   in namespace FTensor. */

#include <cmath>
#include <complex>
#ifdef FTENSOR_DEBUG
#include <iostream>
#endif
#include "Layout.h"
namespace FTensor {
  template <class T> class Tensor0;

  template <class T, int Dim> class Tensor1;
  template<class A, class T, int Dim, char i> class Tensor1_Expr;

  template <class T, int Dim1, int Dim2, Layout layout=column_major>
  class Tensor2;
  template<class A, class T, int Dim1, int Dim2, char i, char j>
  class Tensor2_Expr;
  template<class A, class T, int N>
  class Tensor2_number_rhs_0;
  template<class A, class T, int N>
  class Tensor2_number_rhs_1;

  template <class T, int Dim> class Tensor2_symmetric;
  template<class A, class T, int Dim, char i, char j>
  class Tensor2_symmetric_Expr;

  template <class A, class T, int Dim0, int Dim1, int Dim2,
    char i, char j, char k> class Tensor3_Expr;

  template <class T, int Dim01, int Dim2> class Tensor3_dg;
  template <class A, class T, int Dim01, int Dim2, char i, char j, char k>
  class Tensor3_dg_Expr;
  template<class A, class T, int N> class Tensor3_dg_number_rhs_0;
  template<class A, class T, int N> class Tensor3_dg_number_rhs_2;
  template<class A, class T, int N1, int N2> class Tensor3_dg_number_rhs_01;
  template<class A, class T, int N1, int N2> class Tensor3_dg_number_rhs_12;

  template <class T, int Dim0, int Dim12> class Tensor3_christof;
  template <class A, class T, int Dim0, int Dim12, char i, char j, char k>
  class Tensor3_christof_Expr;

  template <class T, int Dim0, int Dim12> class Tensor3_antisymmetric;
  template <class A, class T, int Dim0, int Dim12, char i, char j, char k>
  class Tensor3_antisymmetric_Expr;

  template <class A, class T, int Dim0, int Dim1, int Dim2, int Dim3,
    char i, char j, char k, char l> class Tensor4_Expr;

  template <class T, int Dim> class Tensor4_Riemann;
  template <class A, class T, int Dim, char i, char j, char k, char l>
  class Tensor4_Riemann_Expr;

  template <class T, int Dim01, int Dim23> class Tensor4_ddg;
  template <class A, class T, int Dim01, int Dim23,
    char i, char j, char k, char l> class Tensor4_ddg_Expr;
  template<class A, class T, int N0, int N1>
  class Tensor4_ddg_number_rhs_01;
  template<class A, class T, int N0>
  class Tensor4_ddg_number_rhs_0;

#include "Index.h"
#include "Number.h"
#include "promote.h"
#include "Tensor0.h"
#include "Tensor1.h"
#include "Tensor2.h"
#include "Tensor2_symmetric.h"
#include "Tensor3/Tensor3_Expr.h"
#include "Tensor3/Tensor3_contracted.h"
#include "Tensor3_dg.h"
#include "Tensor3_christof.h"
#include "Tensor3_antisymmetric.h"
#include "Tensor4/Tensor4_Expr.h"
#include "Tensor4_ddg.h"
#include "Tensor4_Riemann.h"
}
