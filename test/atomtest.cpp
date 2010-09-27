#include "atom.h"

int main () {

	Atom a ("O", VecR(1.0,2.0,3.0));
	Atom b ("H", VecR(2.0,3.0,1.0));

	a.Print();
	b.Print();



	return 0;
}
