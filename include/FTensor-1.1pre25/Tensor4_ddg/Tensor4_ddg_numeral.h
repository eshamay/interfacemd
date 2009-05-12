/* This is for expressions where a number is used for one or more
   slots, and an index for the others, yielding a
   Tensor2_symmetric_Expr or Tensor3_dg_Expr. */

/* First and second slot. */

template<class A, class T>
class Tensor4_ddg_numeral_01
{
  A iterA;
  const int N1, N2;
public:
  T operator()(const int N3, const int N4) const
  {
    return iterA(N1,N2,N3,N4);
  }
  Tensor4_ddg_numeral_01(const A &a, const int NN1, const int NN2)
    :iterA(a), N1(NN1), N2(NN2) {}
};

/* Third and fourth slot. */

template<class A, class T>
class Tensor4_ddg_numeral_23
{
  A iterA;
  const int N3, N4;
public:
  T operator()(const int N1, const int N2) const
  {
    return iterA(N1,N2,N3,N4);
  }
  Tensor4_ddg_numeral_23(const A &a, const int NN3, const int NN4)
    :iterA(a), N3(NN3), N4(NN4) {}
};

/* First and third slot. */

template<class A, class T>
class Tensor4_ddg_numeral_02
{
  A iterA;
  const int N1, N3;
public:
  T operator()(const int N2, const int N4) const
  {
    return iterA(N1,N2,N3,N4);
  }
  Tensor4_ddg_numeral_02(const A &a, const int NN1, const int NN3)
    :iterA(a), N1(NN1), N3(NN3) {}
};

/* Only first slot */

/* I think I have to be careful to make sure that the indices are
   matched correctly, because the expression that calls it is a
   Tensor3_dg. */

template<class A, class T>
class Tensor4_ddg_numeral_0
{
  A iterA;
  const int N1;
public:
  T operator()(const int N2, const int N3, const int N4) const
  {
    return iterA(N1,N4,N2,N3);
  }
  Tensor4_ddg_numeral_0(const A &a, const int NN1):iterA(a), N1(NN1) {}
};

/* Second, third and fourth slot. */

template<class A, class T>
class Tensor4_ddg_numeral_123
{
  A iterA;
  const int N1,N2,N3;
public:
  T operator()(const int N0) const
  {
    return iterA(N0,N1,N2,N3);
  }
  Tensor4_ddg_numeral_123(const A &a, const int NN1, const int NN2,
			  const int NN3) : iterA(a), N1(NN1), N2(NN2), N3(NN3)
  {}
};
