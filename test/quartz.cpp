//#include "../watersystem.h"
#include "../xyzsystem.h"
#include "../pdbfile.h"

using namespace std;

int main () {
	XYZFile xyz ("quartz_unitcell.xyz");
	XYZFile pds ("pds.xyz");
	Atom::Size (VecR(90.0, 90.0, 90.0));

	// create the quartz unitcell
	Molecule uc;
	uc.Name("qtz");
	RUN (xyz) {
		uc.AddAtom(xyz[i]);
	}

	//make copies of the quartz unitcell and shift them around
	VecR x_shift (2.4567,	-4.25512922,	0.0);
	VecR y_shift (2.4567,	 4.25512922,	0.0);
	VecR z_shift (0.0,	 0.0,			5.4052);
	VecR x_solo (2.4567,0.0,0.0);
	VecR y_solo (0.0,4.25512922,0.0);

	vector<Molecule *> mols;
	for (int i = 0; i <6; i++) {
		VecR shift;
		shift += z_shift * (double)i;
	for (int j = 0; j <6; j++) {
		shift += x_shift;
		shift += y_shift;

		// make a copy of the unitcell
		Molecule * copy = new Molecule;
		copy->Name("qtz");
		// make copies of all the atoms in the unitcell
		for (int atom = 0; atom < uc.size(); atom++) {
			Atom * pa = new Atom (*uc[atom]);
			copy->AddAtom(pa);
		}

		// shift the new copy to the next lattice point
		copy->Shift(shift);
		mols.push_back(copy);
	} }

	// create the pds molecule copy
	Molecule pds_copy;
	pds_copy.Name("pds");
	RUN (pds) {
		pds_copy.AddAtom(pds[i]);
	}

	// flip the pds upside down
	VecR x2 (1,0,0);
	VecR y2 (0,1,0);
	VecR z2 (0,0,-1);

	pds_copy.X(x2); pds_copy.Y(y2); pds_copy.Z(z2);
	MatR dcm = pds_copy.DCMToLab(y);

	RUN (pds_copy)
		pds_copy[i]->Position (dcm * pds_copy[i]->Position());

	// find oxygens above which to position a pds
	vector <Atom *> Os;
	RUN (mols) {
		for (int j = 0; j < mols[i]->size(); j++) {
			Atom * o = (*mols[i])[j];
			if (o->Name() != "O") continue;
			if (o->Position()[y] < 3.63) continue;
			Os.push_back(o);
		}
	}

	// place a pds above certain oxygens
	 int iSecret, iGuess;
  	// initialize random seed:
  	srand ( time(NULL) );
	RUN (Os) {

  		// generate random number:
  		iSecret = rand() % 100 + 1;
		if (iSecret > 90) continue;

		// copy the starting molecule
		Molecule * pm = new Molecule;
		pm->Name("pds");
		RUN2 (pds) {
			Atom * pa = new Atom (*pds_copy[j]);
			pm->AddAtom(pa);
		}

		// shift the new pds to where the oxygen is
		Atom * si = pm->GetAtom("Si");
		VecR shift = si->MinVector(Os[i]);
		pm->Shift (shift + VecR(0,2.3,0));
		mols.push_back(pm);
	}

	PDBFile::WritePDB (mols);
return 0;
}
