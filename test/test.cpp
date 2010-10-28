#include "test.h"
int main () {

	Atom O ("O", VecR(0,0,0));
	Atom H1 ("H1", VecR(-.757,0,-.586));
	Atom H2 ("H2", VecR(.757,0,-.586));
	Atom M ("M", VecR(0,0,-.21));

	H1.Position(H1.Position().normalized());
	H2.Position(H2.Position().normalized());

	Water wat;
	wat.AddAtom(&O);
	wat.AddAtom(&H1);
	wat.AddAtom(&H2);
	//wat.AddAtom(&M);
	wat.MolID(0);


	wat.Print();
	wat.SetAtoms();
	wat.UpdateCenterOfMass();
	wat.CenterOfMass().Print();

	std::cout << wat.OH1()->Magnitude() << std::endl;
	std::cout << wat.OH2()->Magnitude() << std::endl;


	return 0;

}

