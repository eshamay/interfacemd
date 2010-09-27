#include "pdbfile.h"

int main () {

	md_files::PDBSystem pdb ("test.pdb");
	for (Mol_it it = pdb.begin_mols(); it != pdb.end_mols(); it++) {
		(*it)->Print();
	}

	return 0;

}
