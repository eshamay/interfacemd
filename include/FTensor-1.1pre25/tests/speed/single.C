#include <iostream.h>
#include <ostream.h>
template<char i>
class Index
{
public:
  Index() {};
};

class Tensor1;

template <class A, char i>
class Tensor1_Expr
{
  A *iter;
public:
  Tensor1_Expr(A *a): iter(a) {}
  double & operator()(const int N)
  {
    return (*iter)(N);
  }
  double operator()(const int N) const
  {
    return (*iter)(N);
  }

  const A operator=(const Tensor1_Expr<Tensor1,'i'> &result)
  {
    cout << "equaling" << endl;
    iter->data0=result(0);
    iter->data1=result(1);
    iter->data2=result(2);
    return *iter;
  }

//    template<class B>
//    const A operator=(const Tensor1_Expr<B,'i'> &result)
//    {
//      cout << "equaling" << endl;
//      iter->data0=result(0);
//      iter->data1=result(1);
//      iter->data2=result(2);
//      return iter;
//    }
};

class Tensor1
{
public:
  double data0, data1, data2;
public:
  Tensor1(double d0, double d1, double d2): data0(d0), data1(d1), data2(d2) {}
  double & operator()(const int N)
  {
    return N==0 ? data0 : (N==1 ? data1 : data2);
  }

  double operator()(const int N) const
  {
    return N==0 ? data0 : (N==1 ? data1 : data2);
  }

  template<char i>
  Tensor1_Expr<Tensor1,i> operator()(const Index<i> index)
  {
    return Tensor1_Expr<Tensor1,i>(this);
  }

  friend ostream& operator<<(ostream& s, const Tensor1 &a);
};

ostream& operator<<(ostream& s, const Tensor1 &a)
{
  return s << a.data0 << " " << a.data1 << " " << a.data2 << " ";
}

int main()
{
  Tensor1 y(0,1,2);
  Tensor1 x(2,3,4);
  const Index<'i'> i;

  y(i)=x(i);

  cout << y << endl;
}
