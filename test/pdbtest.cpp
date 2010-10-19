#include "pdbsystem.h"

int main () {

	md_files::PDBSystem pdb ("test.pdb");
	for (Mol_it it = pdb.begin_mols(); it != pdb.end_mols(); it++) {
		static_cast<Water *>((*it))->H1()->Print();
	}

	return 0;

}
