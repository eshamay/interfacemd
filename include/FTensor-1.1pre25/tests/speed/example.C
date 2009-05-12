#include <iostream>
using namespace std;
#include "../../FTensor.h"
using namespace FTensor;

int main()
{
  Tensor1<double,3> y(0,1,2);
  Tensor1<double,3> x(2,3,4);
  Tensor1<double,3> n(5,6,7);
  const Index<'i',3> i;

  for(int j=0;j<10000000;j++)
    {
      y(i)=x(i)+n(i);
      x(i)=y(i)-n(i);
      n(i)=n(i)+(y(i)-x(i))-(y(i)-x(i))+(y(i)-x(i))-(y(i)-x(i))
	+(y(i)-x(i))-(y(i)-x(i));

      n(i)=n(i)+y(i)-x(i)
 +y(i)-x(i)
 +y(i)-x(i)
 +y(i)-x(i)
 +y(i)-x(i)
 +y(i)-x(i)
 +y(i)-x(i)
 +y(i)-x(i)
 +y(i)-x(i)
 +y(i)-x(i)
 +y(i)-x(i)
 +y(i)-x(i)
 +y(i)-x(i)
 +y(i)-x(i)
 +y(i)-x(i)
 +y(i)-x(i)
 +y(i)-x(i)
 +y(i)-x(i)
 +y(i)-x(i)
 +y(i)-x(i)
 +y(i)-x(i)
 +y(i)-x(i)
 +y(i)-x(i)
 +y(i)-x(i)
 +y(i)-x(i)
 +y(i)-x(i)
+y(i)-x(i)
+y(i)-x(i)
+y(i)-x(i)
+y(i)-x(i)
+y(i)-x(i)
+y(i)-x(i)
+y(i)-x(i)
+y(i)-x(i)
+y(i)-x(i)
+y(i)-x(i)
+y(i)-x(i)
+y(i)-x(i)
+y(i)-x(i)
+y(i)-x(i)
+y(i)-x(i)
+y(i)-x(i)
+y(i)-x(i)
+y(i)-x(i)
+y(i)-x(i)
+y(i)-x(i)
+y(i)-x(i)
+y(i)-x(i)
// +y(i)-x(i)
// +y(i)-x(i)
// +y(i)-x(i)
// +y(i)-x(i)
// +y(i)-x(i)
// +y(i)-x(i)
// +y(i)-x(i)
// +y(i)-x(i)
// +y(i)-x(i)
// +y(i)-x(i)
// +y(i)-x(i)
// +y(i)-x(i)
// +y(i)-x(i)
// +y(i)-x(i)
// +y(i)-x(i)
// +y(i)-x(i)
// +y(i)-x(i)
// +y(i)-x(i)
// +y(i)-x(i)
// +y(i)-x(i)
// +y(i)-x(i)
// +y(i)-x(i)
// +y(i)-x(i)
// +y(i)-x(i)
// +y(i)-x(i)
// +y(i)-x(i)
// +y(i)-x(i)
// 	+y(i)-x(i)
	;

//        n(i)=(y(i)-x(i))*(y(i)-x(i))/(n(i));
// n(i)=(y(i)-x(i))*(y(i)-x(i))*(y(i)-x(i))*(y(i)-x(i))*(y(i)-x(i))*(y(i)-x(i))*(y(i)-x(i))*(y(i)-x(i))*(y(i)-x(i))*(y(i)-x(i))*(y(i)-x(i))*(y(i)-x(i))*(y(i)-x(i))*(y(i)-x(i))*(y(i)-x(i))*(y(i)-x(i))/(n(i)*n(i)*n(i)*n(i)*n(i)*n(i)*n(i)*n(i)*n(i)*n(i)*n(i)*n(i)*n(i)*n(i)*n(i));
    }
  cout << y(0) << " " << y(1) << " " << y(2) << endl;
}
