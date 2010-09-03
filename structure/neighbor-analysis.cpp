#include "neighbor-analysis.h"

namespace md_analysis {


	void so2_closest_O_analyzer::Analysis (system_t& t) {
		t.LoadAll();

		s = Atom::FindByElement (t.sys_atoms, Atom::S);
		closest = t.System()->graph.ClosestAtoms (s, 4, Atom::O, true);

		// every timestep - output the 4 nearest oxygens (including those covalently bound)
		for (bondgraph::distance_vec::const_iterator it = closest.begin(); it != closest.end(); it++) {
			fprintf (t.Output(), "% 12.4f", it->first);
		}
		fprintf(t.Output(),"\n");

		return;
	}


	void so2_closest_H_analyzer::Analysis (system_t& t) {
		t.LoadAll();

		MolPtr mol = Molecule::FindByType(t.sys_mols, Molecule::SO2);
		so2 = new SulfurDioxide(mol);
		so2->SetAtoms();

		// every timestep output the distance to the 3 Hydrogens nearest to both of the SO2 oxygens
		o1 = so2->O1();
		closest = t.System()->graph.ClosestAtoms (o1, 3, Atom::H);
		for (bondgraph::distance_vec::const_iterator it = closest.begin(); it != closest.end(); it++) {
			fprintf (t.Output(), "% 12.4f", it->first);
		}

		o2 = so2->O2();
		closest = t.System()->graph.ClosestAtoms (o2, 3, Atom::H);
		for (bondgraph::distance_vec::const_iterator it = closest.begin(); it != closest.end(); it++) {
			fprintf (t.Output(), "% 12.4f", it->first);
		}
		fprintf(t.Output(), "\n");

		delete so2;

		return;
	}

	void so2_hbond_factor_analyzer::Analysis (system_t& t) {
		t.LoadAll();

		MolPtr mol = Molecule::FindByType(t.sys_mols, Molecule::SO2);
		so2 = new SulfurDioxide(mol);
		so2->SetAtoms();

		o1 = so2->O1();
		o2 = so2->O2();


		// find out the hbond-distance factor for each O atom. 
		// The Distance Factor is a dimensionless number that represents the position of a hydrogen between two oxygens.
		// An H exactly in the middle of 2 oxygens will have Q=0.0
		// an H sitting right on top of the water (target) oxygen will have a value of 1.0, and sitting right on top of the source (SO2) oxygen will have a value of -1.0
		// Q = (dOH1 - dOH2)/D
		// dOH1 = distance of H to SO2 Oxygen
		// dOH2 = distance of H to H2O oxygen
		// D = total distance

		bondgraph::distance_pair h1_pair = t.System()->graph.ClosestAtom (o1, Atom::H);
		bondgraph::distance_pair h2_pair = t.System()->graph.ClosestAtom (o2, Atom::H);

		AtomPtr target_o1 = h1_pair.second->ParentMolecule()->GetAtom(Atom::O);
		AtomPtr target_o2 = h2_pair.second->ParentMolecule()->GetAtom(Atom::O);

		double target_distance1 = t.System()->graph.Distance(target_o1, h1_pair.second);
		double target_distance2 = t.System()->graph.Distance(target_o2, h2_pair.second);

		double Q1 = (h1_pair.first - target_distance1)/(h1_pair.first + target_distance1);
		double Q2 = (h2_pair.first - target_distance2)/(h2_pair.first + target_distance2);

		// data output 
		fprintf (t.Output(), "% 8d % 8.4f % 8.4f\n", t.timestep, Q1, Q2);

		delete so2;

		return;
	}


}		// namespace md_analysis
