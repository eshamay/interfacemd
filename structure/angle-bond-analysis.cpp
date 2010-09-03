#include "angle-bond-analysis.h"

namespace md_analysis {

	void so2_angle_bond_analyzer::Analysis (system_t& t) {
		t.LoadAll();

		MolPtr mol = Molecule::FindByType(t.sys_mols, Molecule::SO2);
		so2 = new SulfurDioxide(mol);
		so2->SetAtoms();

		angle = so2->Angle();
		// output the two S-O bondlengths and the SO2 angle for each timestep
		fprintf (t.Output(), "% 12.4f % 12.4f % 12.4f\n", so2->SO1().norm(), so2->SO2().norm(), acos(angle)*180.0/M_PI);

		delete so2;

		return;
	}

}	// namespace md_analysis
