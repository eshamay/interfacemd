/* This is for expressions where a number is used for one slot, and
   an index for the other, yielding a Tensor2_symmetric_Expr. */

template<class A, class T, int N>
class Tensor3_christof_number_0
{
  const A iterA;
public:
  T operator()(const int N1, const int N2) const
  {
    return iterA(N,N1,N2);
  }
  Tensor3_christof_number_0(const A &a): iterA(a) {}
};

