#ifndef NEIGHBOR_ANALYSIS_H_
#define NEIGHBOR_ANALYSIS_H_

#include "analysis.h"

namespace md_analysis {

	template <typename T>
	class so2_closest_water_map : public AnalysisSet< Analyzer<T> > {
		public:
			typedef Analyzer<T> system_t;
			so2_closest_water_map () :
				AnalysisSet<system_t> (
						std::string("SO2 closest neighbor 3D histogram mapping analysis"),
						std::string ("temp.dat")),
				// create the histograms with distances up to 10.0 angstroms away (in both positive & negative directions)
				xy_histo (std::make_pair(-10.0,-10.0), std::make_pair(10.0,10.0), std::make_pair(0.1,0.1)),
				yz_histo (std::make_pair(-10.0,-10.0), std::make_pair(10.0,10.0), std::make_pair(0.1,0.1)),
				xz_histo (std::make_pair(-10.0,-10.0), std::make_pair(10.0,10.0), std::make_pair(0.1,0.1)) { }

			void Setup (system_t& t) {
				t.LoadAll();
				_so2_mol_id = t.SystemParameterLookup ("analysis.reference-molecule-id");
				MolPtr mol = Molecule::FindByID (t.sys_mols, _so2_mol_id);
				so2 = new SulfurDioxide(mol);
				so2->SetAtoms();
				xy_file = fopen("so2-neighbormap.xy.dat", "w");
				yz_file = fopen("so2-neighbormap.yz.dat", "w");
				xz_file = fopen("so2-neighbormap.xz.dat", "w");
			}

			virtual ~so2_closest_water_map () {
				fclose(xy_file);
				fclose(yz_file);
				fclose(xz_file);
				delete so2;
			}

			void Analysis (system_t& t);
			void DataOutput (system_t& t);

		protected:
			int _so2_mol_id;
			SulfurDioxide * so2;
			// three histograms that hold atom location in the 3 different planes of the molecular_frame
			histogram_utilities::Histogram2D<double>	xy_histo, yz_histo, xz_histo;
			FILE *xy_file, *yz_file, *xz_file;
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
			t.LoadWaters();

			// gather all the system waters
			Water_ptr_vec all_wats, analysis_wats;
			for (Mol_it it = t.int_wats.begin(); it != t.int_wats.end(); it++) {
				WaterPtr wat (new Water(*(*it)));
				wat->SetAtoms();
				all_wats.push_back(wat);
			}
			std::copy (all_wats.begin(), all_wats.end(), std::back_inserter(analysis_wats));

			// find the closest waters by sorting through the water list
			// sort based on the water oxygen distance to the so2-S
			//std::sort(analysis_wats.begin(), analysis_wats.end(), Analyzer<T>::molecule_reference_distance_pred(so2));
			
			// clear all but the Hs in the list
			Atom::KeepByElement (t.int_atoms, Atom::H);
			// Find the Hs that are closest to the so2-O1
			std::sort(t.int_atoms.begin(), t.int_atoms.end(), Analyzer<T>::atomic_reference_distance_pred (so2->O1()));

			// generate the direction cosine matrix from the system to the so2-molecular frame
			so2->SetBisectorAxes();
			so2->DCMToLab();

			// for the closest waters, add their locations into the histograms
			VecR r;
			for (int i = 0; i < 2; i++) {
				// find the distance from the S to the water oxygens
				//r = analysis_wats[i]->ReferencePoint() - so2->ReferencePoint();

				r = t.int_atoms[i]->Position() - so2->ReferencePoint();
				//
				// find out what r looks like in the molecular frame;
				r = so2->DCM().transpose() * r;
				//if (r.norm() < 6.0)
					//fprintf (xy_file, "% .3f % .3f % .3f\n", r[x], r[y], r[z]);
				// now bin the position in each histogram
				xy_histo (r[x], r[y]);
				yz_histo (r[y], r[z]);
				xz_histo (r[x], r[z]);

				/*
				r = analysis_wats[i]->H2()->Position() - so2->ReferencePoint();
				// find out what r looks like in the molecular frame;
				r = so2->DCM().transpose() * r;
				//if (r.norm() < 6.0)
					//fprintf (xy_file, "% .3f % .3f % .3f\n", r[x], r[y], r[z]);

				// now bin the position in each histogram
				xy_histo (r[x], r[y]);
				yz_histo (r[y], r[z]);
				xz_histo (r[x], r[z]);
				*/
			}
		}

	template <typename T>
		void so2_closest_water_map<T>::DataOutput (system_t& t) {
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
		}

}	// namespace md_analysis

#endif
