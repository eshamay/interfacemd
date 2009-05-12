#include "../../FTensor.h"
using namespace FTensor;

inline void func1(Tensor1<double,3> &y, Tensor1<double,3> &a1)
{
  const Index<'i',3> i;

  y(i)+=a1(i);
  a1(i)*=0.1;
  return;
}

inline void func2(Tensor1<double,3> &y, Tensor1<double,3> &a1, Tensor1<double,3> &a2)
{
  const Index<'i',3> i;

  y(i)+=a1(i)
    + 2*a2(i)
    ;
  a1(i)*=0.1;
  a2(i)*=0.2;
  return;
}

inline void func3(Tensor1<double,3> &y, Tensor1<double,3> &a1, Tensor1<double,3> &a2, Tensor1<double,3> &a3)
{
  const Index<'i',3> i;
  const Index<'j',3> j;

  y(i)+=a1(i)
    + 2*a2(i)
    + 3*(a1(j)*a2(j))*a3(i)
    ;
  a1(i)*=0.1;
  a2(i)*=0.2;
  a3(i)*=0.3;
  return;
}

inline void func4(Tensor1<double,3> &y, Tensor1<double,3> &a1, Tensor1<double,3> &a2, Tensor1<double,3> &a3,
		  Tensor1<double,3> &a4)
{
  const Index<'i',3> i;
  const Index<'j',3> j;
  const Index<'k',3> k;

  y(i)+=a1(i)
    + 2*a2(i)
    + 3*(a1(j)*a2(j))*a3(i)
    + 4*(a1(j)*a3(j))*(a2(k)*a2(k))*a4(i)
    ;
  a1(i)*=0.1;
  a2(i)*=0.2;
  a3(i)*=0.3;
  a4(i)*=0.4;
  return;
}

inline void func5(Tensor1<double,3> &y, Tensor1<double,3> &a1, Tensor1<double,3> &a2, Tensor1<double,3> &a3,
		  Tensor1<double,3> &a4, Tensor1<double,3> &a5)
{
  const Index<'i',3> i;
  const Index<'j',3> j;
  const Index<'k',3> k;

  y(i)+=a1(i)
    + 2*a2(i)
    + 3*(a1(j)*a2(j))*a3(i)
    + 4*(a1(j)*a3(j))*(a2(k)*a2(k))*a4(i)
    + 5*(a1(j)*a4(j))*(a2(k)*a3(k))*a5(i)
    ;
  a1(i)*=0.1;
  a2(i)*=0.2;
  a3(i)*=0.3;
  a4(i)*=0.4;
  a5(i)*=0.5;
  return;
}

inline void func6(Tensor1<double,3> &y, Tensor1<double,3> &a1, Tensor1<double,3> &a2, Tensor1<double,3> &a3,
		  Tensor1<double,3> &a4, Tensor1<double,3> &a5, Tensor1<double,3> &a6)
{
  const Index<'i',3> i;
  const Index<'j',3> j;
  const Index<'k',3> k;
  const Index<'l',3> l;

  y(i)+=a1(i)
    + 2*a2(i)
    + 3*(a1(j)*a2(j))*a3(i)
    + 4*(a1(j)*a3(j))*(a2(k)*a2(k))*a4(i)
    + 5*(a1(j)*a4(j))*(a2(k)*a3(k))*a5(i)
    + 6*(a1(j)*a5(j))*(a2(k)*a4(k))*(a3(l)*a3(l))*a6(i)
    ;
  a1(i)*=0.1;
  a2(i)*=0.2;
  a3(i)*=0.3;
  a4(i)*=0.4;
  a5(i)*=0.5;
  a6(i)*=0.6;
  return;
}

inline void func7(Tensor1<double,3> &y, Tensor1<double,3> &a1, Tensor1<double,3> &a2, Tensor1<double,3> &a3,
		  Tensor1<double,3> &a4, Tensor1<double,3> &a5, Tensor1<double,3> &a6, Tensor1<double,3> &a7)
{
  const Index<'i',3> i;
  const Index<'j',3> j;
  const Index<'k',3> k;
  const Index<'l',3> l;

  y(i)+=a1(i)
    + 2*a2(i)
    + 3*(a1(j)*a2(j))*a3(i)
    + 4*(a1(j)*a3(j))*(a2(k)*a2(k))*a4(i)
    + 5*(a1(j)*a4(j))*(a2(k)*a3(k))*a5(i)
    + 6*(a1(j)*a5(j))*(a2(k)*a4(k))*(a3(l)*a3(l))*a6(i)
    + 7*(a1(j)*a6(j))*(a2(k)*a5(k))*(a3(l)*a4(l))*a7(i)
    ;
  a1(i)*=0.1;
  a2(i)*=0.2;
  a3(i)*=0.3;
  a4(i)*=0.4;
  a5(i)*=0.5;
  a6(i)*=0.6;
  a7(i)*=0.7;
  return;
}

inline void func8(Tensor1<double,3> &y, Tensor1<double,3> &a1, Tensor1<double,3> &a2, Tensor1<double,3> &a3,
		  Tensor1<double,3> &a4, Tensor1<double,3> &a5, Tensor1<double,3> &a6, Tensor1<double,3> &a7,
		  Tensor1<double,3> &a8)
{
  const Index<'i',3> i;
  const Index<'j',3> j;
  const Index<'k',3> k;
  const Index<'l',3> l;
  const Index<'m',3> m;

  y(i)+=a1(i)
    + 2*a2(i)
    + 3*(a1(j)*a2(j))*a3(i)
    + 4*(a1(j)*a3(j))*(a2(k)*a2(k))*a4(i)
    + 5*(a1(j)*a4(j))*(a2(k)*a3(k))*a5(i)
    + 6*(a1(j)*a5(j))*(a2(k)*a4(k))*(a3(l)*a3(l))*a6(i)
    + 7*(a1(j)*a6(j))*(a2(k)*a5(k))*(a3(l)*a4(l))*a7(i)
    + 8*(a1(j)*a7(j))*(a2(k)*a6(k))*(a3(l)*a5(l))*(a4(m)*a4(m))*a8(i)
    ;
  a1(i)*=0.1;
  a2(i)*=0.2;
  a3(i)*=0.3;
  a4(i)*=0.4;
  a5(i)*=0.5;
  a6(i)*=0.6;
  a7(i)*=0.7;
  a8(i)*=0.8;
  return;
}

inline void func9(Tensor1<double,3> &y, Tensor1<double,3> &a1, Tensor1<double,3> &a2, Tensor1<double,3> &a3,
		  Tensor1<double,3> &a4, Tensor1<double,3> &a5, Tensor1<double,3> &a6, Tensor1<double,3> &a7,
		  Tensor1<double,3> &a8, Tensor1<double,3> &a9)
{
  const Index<'i',3> i;
  const Index<'j',3> j;
  const Index<'k',3> k;
  const Index<'l',3> l;
  const Index<'m',3> m;

  y(i)+=a1(i)
    + 2*a2(i)
    + 3*(a1(j)*a2(j))*a3(i)
    + 4*(a1(j)*a3(j))*(a2(k)*a2(k))*a4(i)
    + 5*(a1(j)*a4(j))*(a2(k)*a3(k))*a5(i)
    + 6*(a1(j)*a5(j))*(a2(k)*a4(k))*(a3(l)*a3(l))*a6(i)
    + 7*(a1(j)*a6(j))*(a2(k)*a5(k))*(a3(l)*a4(l))*a7(i)
    + 8*(a1(j)*a7(j))*(a2(k)*a6(k))*(a3(l)*a5(l))*(a4(m)*a4(m))*a8(i)
    + 9*(a1(j)*a8(j))*(a2(k)*a7(k))*(a3(l)*a6(l))*(a4(m)*a5(m))*a9(i)
    ;
  a1(i)*=0.1;
  a2(i)*=0.2;
  a3(i)*=0.3;
  a4(i)*=0.4;
  a5(i)*=0.5;
  a6(i)*=0.6;
  a7(i)*=0.7;
  a8(i)*=0.8;
  a9(i)*=0.9;
  return;
}
