#include <iostream>
#include "/common/src/FTensor-1.1pre25/FTensor.h"

using namespace FTensor;
int main () {

	Tensor1<double,3> v (8,5,1);
	const Index<'i',3> i;

	std::cout << v(i) * v(i) << std::endl;

return 0;
}
