#ifndef __ANGLE_BOND_ANALYSIS_H_
#define __ANGLE_BOND_ANALYSIS_H_

#include "analysis.h"
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
						WaterSystem<T>::SystemParameterLookup("analysis.angle-bond-histogram.bondlength-res"),
						WaterSystem<T>::SystemParameterLookup("analysis.angle-bond-histogram.angle-min"), 
						WaterSystem<T>::SystemParameterLookup("analysis.angle-bond-histogram.angle-max"), 
						WaterSystem<T>::SystemParameterLookup("analysis.angle-bond-histogram.angle-res")) { }

			virtual void Analysis (system_t&);

		protected:
			Water * wat;

	};	 // h2o angle & bond histogram analyzer


	class so2_angle_bond_analyzer : public XYZAnalysisSet {
		public:
			virtual ~so2_angle_bond_analyzer () { }
			so2_angle_bond_analyzer () :
				XYZAnalysisSet (
						std::string("[CP2K] SO2 molecular angle and S-O bondlengths"),
						std::string ("so2-angle+bonds.dat")) { }

			void Setup (system_t& t) {
				XYZAnalysisSet::Setup(t);
				fprintf (t.Output(), "so1 so2 theta\n");
			}

			void Analysis (system_t& t);

		protected:
			SulfurDioxide * so2;
			double angle;
	};


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
