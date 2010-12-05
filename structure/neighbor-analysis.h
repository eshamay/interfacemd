#ifndef NEIGHBOR_ANALYSIS_H_
#define NEIGHBOR_ANALYSIS_H_

#include "analysis.h"
#include "so2-system-analysis.h"

namespace md_analysis {

	template <typename T>
	class so2_closest_water_map : public so2_analysis::SO2SystemAnalyzer<T> {
		public:
			typedef typename so2_analysis::SO2SystemAnalyzer<T>::system_t system_t;

			so2_closest_water_map () :
				so2_analysis::SO2SystemAnalyzer<T> (
						std::string("SO2 closest neighbor 3D histogram mapping analysis"),
						std::string ("temp.xyz")) { }

			virtual ~so2_closest_water_map () { }

			void FindSO2 (system_t& t) {
				// find the so2 molecule of interest
				_so2_mol_id = t.SystemParameterLookup ("analysis.reference-molecule-id");
				MolPtr mol = Molecule::FindByID (t.sys_mols, _so2_mol_id);
				this->so2 = new SulfurDioxide(mol);
			}

			void PostSetup (system_t& t) {
				numWatsToProcess = 5;
			}

			void Analysis (system_t& t);
			void DataOutput (system_t& t);

		protected:
			int _so2_mol_id;
			int numWatsToProcess;


			// prints out a single line to the output file comprised of the atomic coordinates of an atom in xyz format
			void PrintXYZLocation (system_t& t, std::string name, VecR r) {
				fprintf (t.Output(), "%6s % .4f % .4f % .4f\n", name.c_str(), r[x], r[y], r[z]);
			}


	};



	class so2_closest_atoms_analyzer : public XYZAnalysisSet {
		public:
			so2_closest_atoms_analyzer (std::string desc, std::string fn) : XYZAnalysisSet (desc,fn) { }
			virtual ~so2_closest_atoms_analyzer () { }
		protected:
			SulfurDioxide * so2;
			AtomPtr s;
			AtomPtr o1,o2;
			bondgraph::distance_vec closest;
	};



	class so2_closest_H_analyzer : public so2_closest_atoms_analyzer {
		public:
			virtual ~so2_closest_H_analyzer () { }
			so2_closest_H_analyzer () :
				so2_closest_atoms_analyzer(
						std::string("[CP2K] SO2 closest hydrogen analysis (reports the distance to the 3 hydrogens closest to each of the SO2 oxygens)"),
						std::string("so2-closest-Hs.dat")) { }

			void Setup (system_t& t) {
				XYZAnalysisSet::Setup(t);
				fprintf (t.Output(), "o11 o12 o13 o21 o22 o23\n");
			}
			void Analysis (system_t& t);
	};



	class so2_closest_O_analyzer : public so2_closest_atoms_analyzer {
		public:
			virtual ~so2_closest_O_analyzer () { }
			so2_closest_O_analyzer () :
				so2_closest_atoms_analyzer(
						std::string ("[CP2K] SO2 closest oxygen analysis (reports the distance to the 4 oxygens closest to the SO2 sulfur)"),
						std::string ("so2-closest-Os.dat")) { }

			void Analysis (system_t& t);
	};


	class so2_closest_OH_analyzer : public so2_closest_atoms_analyzer {
		public:
			virtual ~so2_closest_OH_analyzer () { }
			so2_closest_OH_analyzer () :
				so2_closest_atoms_analyzer(
						std::string ("[CP2K] SO2 neighbor analysis to find the closest oxygens and hydrogens to the SO2. Output is 8 columns - 2x Oxygens closest to S, 3x H's closest to so2-O1, and 3x H's closest to so2-O2"),
						std::string ("so2-closest-Os+Hs.dat")) { }

			void Analysis (system_t& t);
		protected:
			bondgraph::distance_vec closestOs;
			bondgraph::distance_vec closestHs_1;
			bondgraph::distance_vec closestHs_2;
	};


	class so2_hbond_factor_analyzer : public XYZAnalysisSet {
		public:
			virtual ~so2_hbond_factor_analyzer () { }
			so2_hbond_factor_analyzer () :
				XYZAnalysisSet (
						std::string ("[CP2K] H-sharing factor - unitless factor for studying the hydrogen-bond character of neighboring H's"),
						std::string ("so2-hbond-factors.dat")) { }

			void Setup (system_t& t) {
				XYZAnalysisSet::Setup (t);
				fprintf(t.Output(), "timestep q1 q2\n");
			}

			void Analysis (system_t& t);
		protected:
			SulfurDioxide * so2;
			AtomPtr o1,o2;
	};


	template <typename T>
		void so2_closest_water_map<T>::Analysis (system_t& t) {
			// copy back all the atom and water pointers for a new round of analysis
			this->ReloadAnalysisWaters();

			// generate the direction cosine rotation matrix to rotate data from the system to the so2-molecular frame
			this->so2->SetBisectorAxes();
			this->so2->DCMToLab();
			MatR dcm = this->so2->DCM().transpose();


			// This is the reference point for the origin of the mapping
			AtomPtr reference_atom = this->s;
			VecR reference_point = reference_atom->Position();
			// clear all but the atoms of interest in the list
			Atom::KeepByElement (this->analysis_atoms, Atom::O);
			// Find the atoms that are closest to the so2-S
			std::sort(this->analysis_atoms.begin(), this->analysis_atoms.end(), Analyzer<T>::atomic_reference_distance_pred (reference_atom));

			// a list of all the nearest molecules
			std::vector<MolPtr> NearestMols;
			NearestMols.push_back(this->analysis_atoms[0]->ParentMolecule());
			NearestMols.push_back(this->analysis_atoms[1]->ParentMolecule());

			this->ReloadAnalysisWaters();
			MolPtr mol;
			Atom::KeepByElement (this->analysis_atoms, Atom::H);
			std::sort(this->analysis_atoms.begin(), this->analysis_atoms.end(), Analyzer<T>::atomic_reference_distance_pred (this->o1));
			for (int i = 0; i < 2; i++) {
				mol = this->analysis_atoms[i]->ParentMolecule();
				NearestMols.push_back(mol);
			}
			std::sort(this->analysis_atoms.begin(), this->analysis_atoms.end(), Analyzer<T>::atomic_reference_distance_pred (this->o2));
			for (int i = 0; i < 2; i++) {
				mol = this->analysis_atoms[i]->ParentMolecule();
				NearestMols.push_back(mol);
			}

			fprintf (t.Output(), "%d\n\n", NearestMols.size()*3+3);
			VecR r (0.0,0.0,0.0);
			PrintXYZLocation (t, "S", r);
			r = dcm * this->so2->SO1();
			PrintXYZLocation (t, "O", r);
			r = dcm * this->so2->SO2();
			PrintXYZLocation (t, "O", r);

			for (int i = 0; i < NearestMols.size(); i++) {
				mol = NearestMols[i];
				// print out every atom in the water molecule associated with the closest O
				for (Atom_it it = mol->begin(); it != mol->end(); it++) {
					r = MDSystem::Distance (reference_point, (*it)->Position());
					// find out what r looks like in the molecular frame;
					r = dcm * r;
					PrintXYZLocation (t, Atom::Element2String((*it)->Element()), r);
				}
			}
			/*
				 reference_atom = so2->O2();
				 reference_point = reference_atom->Position();
			// Find the atoms that are closest to the so2-S
			std::sort(analysis_atoms.begin(), analysis_atoms.end(), Analyzer<T>::atomic_reference_distance_pred (reference_atom));
			for (int i = 0; i < 2; i++) {
			r = MDSystem::Distance(reference_point, analysis_atoms[i]->Position());
			r = so2->DCM().transpose() * r;
			fprintf (t.Output(), "% .3f % .3f % .3f\n", r[x], r[y], r[z]);
			}
			*/
		}

	template <typename T>
		void so2_closest_water_map<T>::DataOutput (system_t& t) {
			/*
				 rewind(xy_file);
				 rewind(yz_file);
				 rewind(xz_file);

			// for each histogram, output the array of row x column populations
			for (double x_val = xy_histo.min.first; x_val < xy_histo.max.first; x_val += xy_histo.resolution.first) {
			for (double y_val = xy_histo.min.second; y_val < xy_histo.max.second; y_val += xy_histo.resolution.second) {
			fprintf (xy_file, "%.3f ", xy_histo.Population(x_val, y_val));
			}
			fprintf (xy_file, "\n");
			}

			for (double y_val = yz_histo.min.first; y_val < yz_histo.max.first; y_val += yz_histo.resolution.first) {
			for (double z_val = yz_histo.min.second; z_val < yz_histo.max.second; z_val += yz_histo.resolution.second) {
			fprintf (yz_file, "%.3f ", yz_histo.Population(y_val, z_val));
			}
			fprintf (yz_file, "\n");
			}

			for (double x_val = xz_histo.min.first; x_val < xz_histo.max.first; x_val += xz_histo.resolution.first) {
			for (double z_val = xz_histo.min.second; z_val < xz_histo.max.second; z_val += xz_histo.resolution.second) {
			fprintf (xz_file, "%.3f ", xz_histo.Population(x_val, z_val));
			}
			fprintf (xz_file, "\n");
			}
			*/
		}

}	// namespace md_analysis

#endif
