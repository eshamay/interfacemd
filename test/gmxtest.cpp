#include "gmxtest.h"

int main () {

	gromacs::GMXSystem<gromacs::XTCFile> sys ("grofile", "xtcfile");

	MolPtr mol = sys.Molecules(0);

	//for (Mol_it it = sys.begin_mols(); it != sys.end_mols(); it++) {
		//(*it)->Print();
	//}
	//sys.Molecules(0)->Print();
	for (int i = 0; i < 10; i++) {
		//mol->operator[](Atom::O)->Print();
		mol->Print();
		sys.LoadNext();
	}

  return 0;

}
