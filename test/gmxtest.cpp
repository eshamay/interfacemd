#include "gmxtest.h"

int main () {

	GMXSystem sys ("grofile", "trrfile");

	//for (Mol_it it = sys.begin_mols(); it != sys.end_mols(); it++) {
		//(*it)->Print();
	//}
	//sys.Molecules(0)->Print();
	for (int i = 0; i < 10000; i++) {
		sys.LoadNext();
		if ((i % 1000)) continue;
		std::cout << i << std::endl;
		sys.Molecules(0)->Print();
	}

  return 0;

}
