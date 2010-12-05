#include "rotate-unitcell.h"

using namespace std;

int main () {
	md_files::PDBSystem pdb ("beta_cristobalite.unitcell.rotated.pdb");

  // create the unitcell
	MolPtr uc = new Molecule();
	for (Atom_it it = pdb.begin(); it != pdb.end(); it++) {
		uc->AddAtom(*it);
	}
	uc->MolID(1);
	uc->Rename("bc");
	uc->FixAtoms();

	// axes that define the frame of the unitcell
	VecR alpha (1,-1,0);
	VecR beta (1,1,-2);
	VecR gamma (1,1,1);

	// Get the direction cosine matrix that will rotate from the unitcell to the Lab frame
	uc->X(alpha);
	uc->Y(beta);
	uc->Z(gamma);
	MatR dcm = uc->DCMToLab();

	// now rotate the position of every atom in the unitcell to the new frame
	//for (Atom_it it = uc->begin(); it != uc->end(); it++) {
		//(*it)->Position(dcm.transpose() * (*it)->Position());
	//}
	

	MolPtr newmol = new Molecule(*uc);

  //make copies of the quartz unitcell and shift them around
	// the shift vectors are the lab-frame vectors that have been rotated into the new unitcell frame.
  //VecR x_shift (dcm.transpose() * VecR(3.583, 0.0, 0.0));
  //VecR y_shift (dcm.transpose() * VecR(0.0, 3.583, 0.0));
  //VecR z_shift (dcm.transpose() * VecR(0.0, 0.0, 3.583));

  VecR x_shift (5.067127, 0.0, 0.0);
  VecR y_shift (0.0, 8.7765217, 0.0);
  VecR z_shift (0.0, 0.0, 6.205938);

	VecR x_new = x_shift;
	VecR y_new = y_shift;
	VecR z_new = z_shift;

	int x_num = 8;
	int y_num = 8;
	int z_num = 3;

	VecR shift = VecR();
	for (int i = 0; i < z_num; i++) {
		for (int j = 0; j < y_num; j++) {
			for (int k = 0; k < x_num; k++) {
				shift += x_new;

				for (Atom_it jt = uc->begin(); jt != uc->end(); jt++) {
					AtomPtr atom = new Atom(**jt);
					atom->Shift(shift);
					newmol->AddAtom(atom);
				}

			}
			shift += y_new;
			shift -= x_new * x_num;
		}
		shift += z_new;
		shift -= y_new * y_num;
	}

	Mol_ptr_vec mols;
	mols.push_back(newmol);
	md_files::PDBSystem::WritePDB (mols);
	//Atom_ptr_vec atoms = newmol->Atoms();
	//XYZFile::WriteXYZ (atoms);


	return 0;
}
