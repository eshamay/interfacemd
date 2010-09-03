#ifndef DIPOLE_ANALYSIS_H_
#define DIPOLE_ANALYSIS_H_

#include "analysis.h"
#include "histogram-analysis.h"

namespace md_analysis {

	template <typename T>
	class h2o_dipole_magnitude_histogram_analyzer : public histogram_analyzer<T> {
		public:
			typedef typename histogram_analyzer<T>::system_t system_t;

			h2o_dipole_magnitude_histogram_analyzer () :
				histogram_analyzer<T> (
						std::string ("Generate a histogram of H2O dipole moment magnitudes"),
						std::string ("h2o-dipole-magnitude-histogram.dat")) { }

			void Analysis (system_t& t);
	};	// H2O dipole magnitude histogram


	template <typename T>
	void h2o_dipole_magnitude_histogram_analyzer<T>::Analysis (system_t& t) {
		t.LoadWaters();

		std::for_each (t.int_wats.begin(), t.int_wats.end(), MDSystem::CalcWannierDipole);

		for (Mol_it it = t.int_wats.begin(); it != t.int_wats.end(); it++) {
			this->values.push_back((*it)->Dipole().Magnitude());
		}

		return;
	}

}	// namespace md_analysis

#endif
