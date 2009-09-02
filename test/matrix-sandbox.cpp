#include "../matrixr.h"
#include "../ambersystem.h"

int main () {

	AmberSystem sys ("prmtop", "mdcrd", "mdvel");

	Water * wat;
	RUN (sys.Molecules()) {
		if (sys.Molecules(i)->Name() != "h2o") continue;
		wat = static_cast<Water *>(sys.Molecules(i));
		break;
	}

	wat->Print();
	wat->SetAtoms();

	MatR D = wat->DCMToLabMorita(z,1);
	MatR Dt = D.Transpose();
	printf ("the D\n");
	D.Print();

	MatR dcm2 = wat->DCMToLabMorita(z,2);
	printf ("the dcm2\n");
	dcm2.Print();

	VecR oh1 = *(wat->OH1());
	printf ("the oh1\n");
	oh1.Print();

	VecR oh2 = *(wat->OH2());
	printf ("the oh2\n");
	oh2.Print();

	printf ("oh1 in the oh1 frame\n");
	(Dt*oh1).Print();

	printf ("oh2 in the oh1 frame\n");
	(Dt*oh2).Print();

	printf ("oh2 in the oh2 frame\n");
	(dcm2.Transpose()*oh2).Print();

	MatR m;
	m.Set(0,0,1.0);
	m.Set(1,1,2.0);
	m.Set(2,2,3.0);
	VecR Z (0,0,1);


	printf ("m*Z\n");
	(m*Z).Print();

	MatR ml = D*m*Dt;
	VecR Zl = D*Z;
	(Dt*ml*Zl).Print();


return 0;
}
