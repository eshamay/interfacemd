#include "so2-structure-analysis.h"

void StructureAnalyzer::PromptForAnalysisFunction () {

  analyses.push_back (new h2o_angle_bond_histogram_analyzer);
  analyses.push_back (new so2_angle_bond_histogram_analyzer);
  analyses.push_back (new so2_angle_bond_analyzer);
  analyses.push_back (new so2_closest_H_analyzer);
  analyses.push_back (new so2_closest_O_analyzer);
  analyses.push_back (new so2_hbond_factor_analyzer);
  analyses.push_back (new h2o_dipole_magnitude_histogram_analyzer);


  printf ("Choose the system analysis to perform from the list below\n\n");
  int choice = 0;
  for (std::vector<XYZAnalysisSet *>::const_iterator it = analyses.begin(); it != analyses.end(); it++) {
	printf ("\t%d) %s\n", choice+1, (*it)->Description().c_str());
	++choice;
  }
  printf ("analysis choice:  ");
  std::cin >> choice;
  //printf ("\n\nperforming analysis (%d) using output filename \"%s\"\n", choice, analyses[choice-1]->Filename().c_str());

  this->SystemAnalysis(*analyses[choice-1]);

  return;
}



// *************************************************** analyzers ***************************************************
//
//


void histogram_analyzer::DataOutput (system_t& t) {

  rewind(t.Output());

  histogram::histogram_t first_histo = histogram::Histogram (values.begin(), values.end(), 100);
  int first_histo_max = histogram::MaxPopulation (first_histo.begin(), first_histo.end());

  for (unsigned i = 0; i < first_histo.size(); i++) {
		fprintf (t.Output(), "% 8.4f % 8.4f\n",
				first_histo[i].first, (double)first_histo[i].second/first_histo_max);
  }
  return;
}


void double_histogram_analyzer::DataOutput (system_t& t) {

  rewind(t.Output());

  histogram::histogram_t first_histo = histogram::Histogram (values.begin(), values.end(), 200);
  int first_histo_max = histogram::MaxPopulation (first_histo.begin(), first_histo.end());

  histogram::histogram_t second_histo = histogram::Histogram (second_values.begin(), second_values.end(), 200);
  int second_histo_max = histogram::MaxPopulation (second_histo.begin(), second_histo.end());

  for (unsigned i = 0; i < first_histo.size(); i++) {
	fprintf (t.Output(), "% 8.4f % 8.4f % 8.4f % 8.4f\n", 
		first_histo[i].first, (double)first_histo[i].second/first_histo_max, 
		second_histo[i].first, (double)second_histo[i].second/second_histo_max);
  }
  return;
}


void h2o_angle_bond_histogram_analyzer::Analysis (system_t& t) {
  t.LoadWaters();

  for (Mol_it it = t.int_wats.begin(); it != t.int_wats.end(); it++) {
	wat = new Water(*it);
	wat->SetAtoms();

	values.push_back(wat->OH1()->norm());
	values.push_back(wat->OH2()->norm());
	double angle = wat->Angle();
	second_values.push_back(acos(angle)*180.0/M_PI);
	delete wat;
  }

  return;
}


void so2_angle_bond_histogram_analyzer::Analysis (system_t& t) {
  t.LoadAll();

  MolPtr mol = Molecule::FindByType(t.sys_mols, Molecule::SO2);
  so2 = new SulfurDioxide(mol);
  so2->SetAtoms();

  values.push_back(so2->SO1().norm());
  values.push_back(so2->SO2().norm());
  double angle = so2->Angle();
  second_values.push_back(acos(angle)*180.0/M_PI);
  delete so2;

  return;
}


void h2o_dipole_magnitude_histogram_analyzer::Analysis (system_t& t) {
	t.LoadWaters();

	std::for_each (t.int_wats.begin(), t.int_wats.end(), MDSystem::CalcWannierDipole);

	for (Mol_it it = t.int_wats.begin(); it != t.int_wats.end(); it++) {
		values.push_back((*it)->Dipole().Magnitude());
	}

	return;
}



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


int main () {
	StructureAnalyzer sa;

	return 0;
}


