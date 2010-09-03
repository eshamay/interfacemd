#ifndef NEIGHBOR_ANALYSIS_H_
#define NEIGHBOR_ANALYSIS_H_

#include "analysis.h"

namespace md_analysis {


	class so2_closest_atoms_analyzer : public XYZAnalysisSet {
		public:
			so2_closest_atoms_analyzer (std::string desc, std::string fn) : XYZAnalysisSet (desc,fn) { }
			virtual ~so2_closest_atoms_analyzer () { }
		protected:
			SulfurDioxide * so2;
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
		protected:
			AtomPtr o1,o2;
	};



	class so2_closest_O_analyzer : public so2_closest_atoms_analyzer {
		public:
			virtual ~so2_closest_O_analyzer () { }
			so2_closest_O_analyzer () :
				so2_closest_atoms_analyzer(
						std::string ("[CP2K] SO2 closest oxygen analysis (reports the distance to the 4 oxygens closest to the SO2 sulfur)"),
						std::string ("so2-closest-Os.dat")) { }

			void Analysis (system_t& t);
		protected:
			AtomPtr s;
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

}	// namespace md_analysis

#endif
