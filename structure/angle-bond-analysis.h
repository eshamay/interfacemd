#ifndef __ANGLE_BOND_ANALYSIS_H_
#define __ANGLE_BOND_ANALYSIS_H_

#include "histogram-analysis.h"
#include "so2-system-analysis.h"

namespace md_analysis {


	template <typename T>
	class h2o_angle_bond_histogram_analyzer : public double_histogram_analysis<T> {
		public:
			typedef typename double_histogram_analysis<T>::system_t system_t;

			virtual ~h2o_angle_bond_histogram_analyzer() { }

			h2o_angle_bond_histogram_analyzer () :
				double_histogram_analysis <T> (
						std::string ("H2O H-O-H angle and O-H bondlength histograms"),
						std::string ("h2o-bond-angle-histograms.dat"),
						WaterSystem<T>::SystemParameterLookup("analysis.angle-bond-histogram.bondlength-min"), 
						WaterSystem<T>::SystemParameterLookup("analysis.angle-bond-histogram.bondlength-max"), 
						WaterSystem<T>::SystemParameterLookup("analysis.angle-bond-histogram.angle-min"), 
						WaterSystem<T>::SystemParameterLookup("analysis.angle-bond-histogram.angle-max"),
						WaterSystem<T>::SystemParameterLookup("analysis.angle-bond-histogram.number-of-bins")) { }

			virtual void Analysis (system_t&);

		protected:
			Water * wat;

	};	 // h2o angle & bond histogram analyzer


	template <typename T>
	class so2_angle_bond_analyzer : public so2_analysis::SO2SystemAnalyzer<T> {
		public:
			typedef typename so2_analysis::SO2SystemAnalyzer<T>::system_t system_t;

			so2_angle_bond_analyzer () :
				so2_analysis::SO2SystemAnalyzer<T> (
						std::string("SO2 molecular angle and S-O bondlengths"),
						std::string ("so2-angle+bonds.dat")) { }

			virtual ~so2_angle_bond_analyzer () { }

			void PostSetup (system_t& t) {
				//AnalysisSet<system_t>::Setup(t);
				// Ox-H-n == the distance from a given so2-O to the nth closest h2o-H
				// S-O-n == distance from so2-S to the nth closest h2o-O
				fprintf (t.Output(), "step location distance SO-1 SO-2 OSO-theta Bisector-theta S-O-1 S-O-2 S-O-3 O1-H-1 O1-H-2 O1-H-3 O2-H-1 O2-H-2 O2-H-3\n");

				// grab the first location of the so2 as a reference for removing all the extraneous waters above it in later analyses
				starting_point = system_t::Position(this->s);
			}

			void Analysis (system_t& t);

		protected:
			double angle;
			double starting_point;	// the starting location of the so2 (along the reference axis)
	};




	template <typename T>
		void so2_angle_bond_analyzer<T>::Analysis (system_t& t) {

			// calculate the OSO angle
			double oso_angle = acos(this->so2->Angle())*180.0/M_PI;
			VecR bisector = this->so2->Bisector();
			VecR Y = Vector3d::UnitY();
			// and the angle of the bisector to the system normal
			double system_angle = acos(bisector < Y)*180.0/M_PI;

			this->ReloadAnalysisWaters ();

			// get rid of everything above the so2
			this->analysis_wats.erase(
					remove_if(this->analysis_wats.begin(), this->analysis_wats.end(), Molecule::MoleculeAbovePosition(starting_point, t.axis)), this->analysis_wats.end());

			// sort the waters by position along the reference axis - first waters are lowest, last are highest
			std::sort (this->analysis_wats.begin(), this->analysis_wats.end(), system_t::molecule_position_pred(Atom::O));
			int numWats = 20;				// number of waters to use for calculating the location of the "top" of the water surface
			double location = 0.0;
			for (Wat_it it = this->analysis_wats.end() - 1; it != this->analysis_wats.end() - numWats; it--) {
				location += system_t::Position((*it)->O());
			}
			location /= numWats;

			// distance from the top of the water surface to the sulfur of the so2
			double distance_to_location = system_t::Position(this->s) - location;

			// sort the waters by distance: h2o-O to the so2-S - the first waters in the vector will be closest to the SO2
			Atom::KeepByElement(t.int_atoms, Atom::O);

			std::sort(t.int_atoms.begin(), t.int_atoms.end(), 
					Analyzer<T>::atomic_reference_distance_pred(this->so2->S()));

			// grab the distances from so2-S to closest h2o-Os
			double so_1 = MDSystem::Distance(this->so2->S(), t.int_atoms[0]).norm();
			double so_2 = MDSystem::Distance(this->so2->S(), t.int_atoms[1]).norm();
			double so_3 = MDSystem::Distance(this->so2->S(), t.int_atoms[2]).norm();

			// sort the water atoms by distance h2o-H to so2-O1
			t.LoadWaters();
			Atom::KeepByElement(t.int_atoms, Atom::H);
			std::sort(t.int_atoms.begin(), t.int_atoms.end(), Analyzer<T>::atomic_reference_distance_pred(this->so2->GetAtom("O1")));

			// grab the distances from so2-S to closest h2o-Os
			double oh1_1 = MDSystem::Distance(this->so2->O1(), t.int_atoms[0]).norm();
			double oh1_2 = MDSystem::Distance(this->so2->O1(), t.int_atoms[1]).norm();
			double oh1_3 = MDSystem::Distance(this->so2->O1(), t.int_atoms[2]).norm();

			// sort the water atoms by distance h2o-H to so2-O1
			t.LoadWaters();
			Atom::KeepByElement(t.int_atoms, Atom::H);
			std::sort(t.int_atoms.begin(), t.int_atoms.end(), Analyzer<T>::atomic_reference_distance_pred(this->so2->GetAtom("O2")));

			// grab the distances from so2-S to closest h2o-Os
			double oh2_1 = MDSystem::Distance(this->so2->O2(), t.int_atoms[0]).norm();
			double oh2_2 = MDSystem::Distance(this->so2->O2(), t.int_atoms[1]).norm();
			double oh2_3 = MDSystem::Distance(this->so2->O2(), t.int_atoms[2]).norm();

			// find the fixed oxygen that is tethered for the steered MD
			//AtomPtr fixed_o = Atom::FindByID(t.int_atoms, 1161);

			// output the distance and the two S-O bondlengths and the SO2 oso_angle for each timestep
			fprintf (t.Output(), "%d % 12.4f % 12.4f % 12.4f % 12.4f % 12.4f % 12.4f % 12.4f % 12.4f % 12.4f % 12.4f % 12.4f % 12.4f % 12.4f % 12.4f % 12.4f\n", 
					t.Timestep(),
					location,
					distance_to_location,
					this->so2->SO1().norm(), this->so2->SO2().norm(), 
					oso_angle, 
					system_angle,
					so_1, so_2, so_3,
					oh1_1, oh1_2, oh1_3,
					oh2_1, oh2_2, oh2_3);

			return;
		}


	template <typename T>
		class so2_angle_bond_histogram_analyzer : public double_histogram_analyzer<T> {
			public:
				typedef typename double_histogram_analyzer<T>::system_t system_t;

				virtual ~so2_angle_bond_histogram_analyzer() { }
				so2_angle_bond_histogram_analyzer() :
					double_histogram_analyzer<T> (
							std::string("SO2 O-S-O angle and S-O bondlength histograms"),
							std::string("so2-bond-angle-histograms.dat")) { }

				void Analysis (system_t& t);
			protected:
				SulfurDioxide * so2;
		};



	template <typename T>
		void h2o_angle_bond_histogram_analyzer<T>::Analysis (system_t& t) {
			t.LoadWaters();

			for (Mol_it it = t.int_wats.begin(); it != t.int_wats.end(); it++) {
				wat = new Water(*it);
				wat->SetAtoms();

				this->histogram(wat->OH1()->Magnitude());
				this->histogram(wat->OH2()->Magnitude());

				// calculate the angle of the water molecule
				double angle = wat->Angle();
				this->histogram_2(acos(angle)*180.0/M_PI);

				delete wat;
			}

			return;
		}



	template <typename T>
		void so2_angle_bond_histogram_analyzer<T>::Analysis (system_t& t) {
			t.LoadAll();

			MolPtr mol = Molecule::FindByType(t.sys_mols, Molecule::SO2);
			so2 = new SulfurDioxide(mol);
			so2->SetAtoms();

			this->values.push_back(so2->SO1().norm());
			this->values.push_back(so2->SO2().norm());
			double angle = so2->Angle();
			this->second_values.push_back(acos(angle)*180.0/M_PI);
			delete so2;

			return;
		}


}	// namespace md_analysis

#endif
