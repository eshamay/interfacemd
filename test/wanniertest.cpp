#include "../wannier.h"

int main () {

	WannierFile wan ("wanniers");

	for (int i=0;i<5;i++) {
		wan[i].Print();
	}

return 0;
}
