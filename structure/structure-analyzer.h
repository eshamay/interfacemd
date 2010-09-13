#ifndef TEST_H_
#define TEST_H_

#define EIGEN_MATRIXBASE_PLUGIN "/home/eshamay/md/src/EigenMatrixAddon.h"
#include <Eigen/Core>
USING_PART_OF_NAMESPACE_EIGEN

#include "../analysis.h"
#include "../utility.h"
#include "angle-bond-analysis.h"
#include "dipole-analysis.h"
#include "neighbor-analysis.h"
#include "atomic-density-analysis.h"
#include "sfg/amber-morita2002.h"


typedef std::vector<double> double_vec;
typedef double_vec::const_iterator double_it;

namespace md_analysis {

	template <typename T>
		class StructureAnalyzer : public Analyzer<T>
	{
		public:
			typedef Analyzer<T> system_t;

			StructureAnalyzer () 
				: system_t () 
			{
				LoadSystemAnalyses ();
				PromptForAnalysisFunction(); 
			}

		protected:

			// output a bit of text to stdout and ask the user for a choice as to which type of analysis to perform - then do it.
			void PromptForAnalysisFunction ();
			//! Loads all the analyses that are compatible with the given system-type
			void LoadSystemAnalyses ();
			// The set of possible analyses to perform on a given system
			typedef std::vector<AnalysisSet<system_t> *>	analysis_vec;
			analysis_vec analyses;
	};


	//! Loads all the system analyses that can be performed on Amber systems
	template <>
		void StructureAnalyzer<AmberSystem>::LoadSystemAnalyses () {
			typedef AmberSystem T;
			AnalysisSet<system_t> * a;

			a = new h2o_angle_bond_histogram_analyzer<T>();
			analyses.push_back(a);
			a = new so2_angle_bond_histogram_analyzer<T>();
			analyses.push_back(a);
			a = new h2o_dipole_magnitude_histogram_analyzer<T>();
			analyses.push_back(a);
			a = new atomic_density_analysis<T>();
			analyses.push_back(a);
			a = new morita::AmberMorita2008Analysis();
			analyses.push_back(a);
		}

	template <>
		void StructureAnalyzer<XYZSystem>::LoadSystemAnalyses () {
			typedef XYZSystem T;
			AnalysisSet<system_t> * a;

			a = new h2o_angle_bond_histogram_analyzer<T>();
			analyses.push_back(a);
			a = new so2_angle_bond_histogram_analyzer<T>();
			analyses.push_back(a);
			a = new h2o_dipole_magnitude_histogram_analyzer<T>();
			analyses.push_back(a);
			a = new atomic_density_analysis<T>();
			analyses.push_back(a);
			a = new so2_angle_bond_analyzer();
			analyses.push_back(a);
			a = new so2_closest_H_analyzer();
			analyses.push_back(a);
			a = new so2_closest_O_analyzer();
			analyses.push_back(a);
			a = new so2_hbond_factor_analyzer();
			analyses.push_back(a);
		}


	template <typename T>
		void StructureAnalyzer<T>::PromptForAnalysisFunction () {

			printf ("Choose the system analysis to perform from the list below\n\n");
			int choice = 0;
			for (typename analysis_vec::iterator it = analyses.begin(); it != analyses.end(); it++) {
				printf ("\t%d) %s\n", choice+1, (*it)->Description().c_str());
				++choice;
			}
			printf ("analysis choice:  ");
			std::cin >> choice;
			//printf ("\n\nperforming analysis (%d) using output filename \"%s\"\n", choice, analyses[choice-1]->Filename().c_str());

			this->SystemAnalysis(*analyses[choice-1]);

			return;
		}

}

#endif

