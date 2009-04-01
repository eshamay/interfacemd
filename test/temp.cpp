#include <iostream>
#include "../ambersystem.h"

using namespace std;

int main () {

	AmberSystem sys ("prmtop", "mdcrd", "");

	Molecule * mol = sys.Molecules(456);
	mol->Print();
	cout << endl;
	Water * wat = static_cast<Water *>(mol);
	wat->SetAtoms();
	const VecR * oh1 = wat->OH1();
	const VecR * oh2 = wat->OH2();


	cout << "oh1 in lab frame" << endl;
	oh1->Print();
	cout << endl;
	cout << "oh2 in lab frame" << endl;
	oh2->Print();
	cout << endl;

	VecR dmu_1_1 (-0.058, 0.0, 0.157);
	VecR dmu_2_1 (0.1287, 0.0, -1.070);

	MatR dcm_1_l = wat->DCMToLabMorita(z);
	VecR oh2_1 (0.9428, 0.0, -0.3333);
	VecR oh1_1 (0.0, 0.0, 1.0);

	// dcm to go from the oh2 frame to the oh1 frame
	double dd[9] = {	0.3333,  0.0, 0.9428,
						0.0,	-1.0, 0.0,
						0.9428,	 0.0, -0.3333 
				   };

	MatR dcm_2_1 (dd);

	double dalpha[9] = {	1.539, 0.0, -0.163,
							0.0, 1.656, 0.0,
							-0.163, 0.0, 7.2 };
	MatR alpha (dalpha);

	(dcm_2_1.Transpose() * alpha * dcm_2_1).Print();

return 0;
}
