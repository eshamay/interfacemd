#ifndef AMBER_MORITA2002_H_
#define AMBER_MORITA2002_H_

#include "morita2002.h"


namespace morita {

	class AmberMorita2008Analysis : public Morita2008LookupAnalysis<AmberSystem> {

		public:
			AmberMorita2008Analysis () :
				Morita2002LookupAnalysis<AmberSystem>(
						std::string("Analysis of the waters in an MD system using the technique of Morita/Hynes (2008) with a system produced by the Amber package (or any compliant dataset (prmtop, mdcrd, mdvel)"),
						std::string("amber-morita2008.dat")) { }

		protected:
			//! Selects the waters to be analyzed by loading the entire set of molecules above a particular cutoff location
			void SelectAnalysisWaters ();
			//! sets the dipole moments of all the system waters by means of the parameterized model of Morita/Hynes 2002

	};	// class Amber sfg analyzer

} // namespace morita
#endif
