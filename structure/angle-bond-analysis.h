#ifndef __ANGLE_BOND_ANALYSIS_H_
#define __ANGLE_BOND_ANALYSIS_H_

#include "histogram-analysis.h"

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
	class so2_angle_bond_analyzer : public AnalysisSet< Analyzer<T> > {
		public:
			typedef Analyzer<T> system_t;
			virtual ~so2_angle_bond_analyzer () { }
			so2_angle_bond_analyzer () :
				AnalysisSet<system_t> (
						std::string("SO2 molecular angle and S-O bondlengths"),
						std::string ("so2-angle+bonds.dat")) { }

			//void Setup (system_t& t) {
				//AnalysisSet<system_t>::Setup(t);
				//fprintf (t.Output(), "distance so1 so2 theta\n");
			//}

			void Analysis (system_t& t);

		protected:
			SulfurDioxide * so2;
			double angle;
	};

	template <typename T>
	void so2_angle_bond_analyzer<T>::Analysis (system_t& t) {
		t.LoadWaters();

		// find the so2 in the system
		MolPtr mol = Molecule::FindByType(t.sys_mols, Molecule::SO2);
		so2 = new SulfurDioxide(mol);
		so2->SetAtoms();
		AtomPtr s = so2->S();
		// calculate the OSO angle
		double oso_angle = acos(so2->Angle())*180.0/M_PI;
		VecR bisector = so2->Bisector();
		VecR Y = Vector3d::UnitY();
		// and the angle of the bisector to the system normal
		double system_angle = acos(bisector < Y)*180.0/M_PI;
		
		// gather all the system waters
		Water_ptr_vec all_wats, analysis_wats;
		for (Mol_it it = t.int_wats.begin(); it != t.int_wats.end(); it++) {
			WaterPtr wat (new Water(*(*it)));
			wat->SetAtoms();
			all_wats.push_back(wat);
		}
		std::copy (all_wats.begin(), all_wats.end(), std::back_inserter(analysis_wats));
		// now all_wats has... all the waters, and analysis wats is used to perform some analysis
		
		// sort the waters by distance to the so2
		std::sort(analysis_wats.begin(), analysis_wats.end(), Analyzer<T>::molecule_distance_pred(so2));

		// and keep only the closest handful of waters
		analysis_wats.erase(analysis_wats.begin(), analysis_wats.end() - 100);

		for (Wat_it it = analysis_wats.begin(); it != analysis_wats.end(); it++) {
			printf ("% 12.4f\n", (so2->ReferencePoint() - (*it)->ReferencePoint()).norm());
		}
		printf ("*****************\n");
		for (Wat_it it = all_wats.begin(); it != all_wats.end(); it++) {
			delete *it;
		}

		/*
		// find the fixed oxygen that is tethered for the steered MD
		AtomPtr fixed_o = Atom::FindByID(t.int_atoms, 1161);

		double distance_to_fixed = (s->Position() - fixed_o->Position()).norm();

		// find the topmost water of the slab
		double pos = 0.0;
		for (Mol_it it = t.sys_mols.begin(); it != t.sys_mols.end(); it++) {
			if ((*it)->Name() != "h2o")
				continue;

			AtomPtr o = (*it)->GetAtom(Atom::O);
			if (o->Position()[y] > pos) {
				pos = o->Position()[y];
			}
		}
		double distance_to_top = s->Position()[y] - pos;

		// output the distance and the two S-O bondlengths and the SO2 oso_angle for each timestep
		//fprintf (t.Output(), "% 12.4f % 12.4f % 12.4f % 12.4f % 12.4f % 12.4f\n", distance_to_fixed, distance_to_top, so2->SO1().norm(), so2->SO2().norm(), oso_angle, system_angle);

		*/
		delete so2;

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
