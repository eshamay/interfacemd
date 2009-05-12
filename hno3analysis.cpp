#include "hno3analysis.h"

//the formula for the contact-ion pair "quality" function is given in the Hynes paper. You can extract it below.
double q_pt (XYZSystem& sys, Molecule * hno3) {

	// we begin by looking establishing the distances between the donor and acceptor oxygens to the proton.
	// We have to first find the acceptor water's oxygen atom

	// here us what we know from the nitric acid
	NitricAcid * na = static_cast<NitricAcid *>(hno3);
	na->SetAtoms();
	Atom * Od = na->GetOH();
	Atom * h = na->GetH();

	// now perform a convoluted routine to grab the oxygen closest to the nitric acid hydrogen on a water molecule
	vector<Atom *> Os (sys.AdjacentAtoms (h, "O"));

	Atom * Oa = (Atom *)NULL;
	RUN (Os) {
		if (Os[i]->Residue() != "h2o" or Os[i]->Residue() != "h3o" or Os[i]->Residue() != "oh") continue;
		Oa = Os[i];
		break;
	}

	if (Oa == (Atom *)NULL) {
		return (-2.0);
		/*
		printf ("while analyzing:\n");
		na->Print();
		printf ("q_pt () :: there was no H-bonded acceptor water found. Now exiting\n");
		exit (1);
		*/
	}

	double distance = *Od - *Oa;
	//printf ("% 10.5f\n", distance);

	// Now let's calculate some values we'll need for the quality function
	double R_oh_a = *Oa - *h;
	double R_oh_d = *Od - *h;

	double d_R_oh_a = R_oh_a - R_OH_H3O_mean;
	double d_R_oh_d = R_oh_d - R_OH_HNO3_mean;

	return ( (d_R_oh_d - d_R_oh_a)/(d_R_oh_d + d_R_oh_a) );
}
