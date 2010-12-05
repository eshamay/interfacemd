#include "quartz-slab.h"

using namespace std;

int main () {
	md_files::PDBSystem pdb ("beta_cristobalite.unitcell.pdb");

  // create the quartz unitcell
	MolPtr uc = new Molecule();
	for (Atom_it it = pdb.begin(); it != pdb.end(); it++) {
		uc->AddAtom(*it);
	}
	uc->MolID(1);
	uc->Rename("bc");
	uc->FixAtoms();

	MolPtr newmol = new Molecule(*uc);

  //make copies of the quartz unitcell and shift them around
  VecR x_shift (3.583, 0.0, 0.0);
  VecR y_shift (0.0, 3.583, 0.0);
	VecR z_shift (0.0, 0.0, 3.583);

	int x_num = 20;
	int y_num = 20;
	int z_num = 20;

	VecR shift = VecR();
	for (int i = 0; i < z_num; i++) {
		for (int j = 0; j < y_num; j++) {
			for (int k = 0; k < x_num; k++) {
				shift += x_shift;

				for (Atom_it jt = uc->begin(); jt != uc->end(); jt++) {
					AtomPtr atom = new Atom(**jt);
					atom->Shift(shift);
					newmol->AddAtom(atom);
				}

			}
			shift += y_shift;
			shift -= x_shift * x_num;
		}
		shift += z_shift;
		shift -= y_shift * y_num;
	}

	Mol_ptr_vec mols;
	mols.push_back(newmol);
	md_files::PDBSystem::WritePDB (mols);
	//Atom_ptr_vec atoms = newmol->Atoms();
	//XYZFile::WriteXYZ (atoms);


	return 0;
}
