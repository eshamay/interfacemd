#include <iostream>
#include "../vecr.h"

using namespace std;

int main () {

	VecR u (1,1,0);
	VecR v (5,3.5,0);
	VecR size (6,4,0);

	VecR w = u.MinVector(v, size);
	w.Print();
	cout << v.MinDistance(u,size);


	

return 0;
}
