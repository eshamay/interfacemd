/* A helper class that automatically returns the internal contraction
   of a rank 3 tensor. */

/* Contracts second and third indices. */

template<class A, class T, int Dim>
class Tensor3_contracted_12
{
  const A iterA;

  template<int Current_Dim>
  T eval(const int N, const Number<Current_Dim> &ND) const
  {
    return iterA(N,Current_Dim-1,Current_Dim-1)
      + eval(N,Number<Current_Dim-1>());
  }
  T eval(const int N, const Number<1> &ND) const
  {
    return iterA(N,0,0);
  }
public:
  T operator()(const int N) const
  {
    return eval(N,Number<Dim>());
  }
  Tensor3_contracted_12(const A &a): iterA(a) {}
};

/* Contracts first and third indices. */

template<class A, class T, int Dim>
class Tensor3_contracted_02
{
  const A iterA;

  template<int Current_Dim>
  T eval(const int N, const Number<Current_Dim> &ND) const
  {
    return iterA(Current_Dim-1,N,Current_Dim-1)
      + eval(N,Number<Current_Dim-1>());
  }
  T eval(const int N, const Number<1> &ND) const
  {
    return iterA(0,N,0);
  }
public:
  T operator()(const int N) const
  {
    return eval(N,Number<Dim>());
  }
  Tensor3_contracted_02(const A &a): iterA(a) {}
};

/* Contracts first and second indices. */

template<class A, class T, int Dim>
class Tensor3_contracted_01
{
  const A iterA;

  template<int Current_Dim>
  T eval(const int N, const Number<Current_Dim> &ND) const
  {
    return iterA(Current_Dim-1,Current_Dim-1,N)
      + eval(N,Number<Current_Dim-1>());
  }
  T eval(const int N, const Number<1> &ND) const
  {
    return iterA(0,0,N);
  }
public:
  T operator()(const int N) const
  {
    return eval(N,Number<Dim>());
  }
  Tensor3_contracted_01(const A &a): iterA(a) {}
};
