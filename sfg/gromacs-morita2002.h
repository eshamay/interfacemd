#ifndef AMBER_MORITA2002_H_
#define AMBER_MORITA2002_H_

#include "morita2002.h"

namespace morita {

	class GMXMorita2008Analysis : public Morita2008LookupAnalysis<GMXSystem> {

		public:
			GMXMorita2008Analysis () :
				Morita2008LookupAnalysis<GMXSystem>(
						std::string("Analysis of the waters in an MD system using the technique of Morita/Hynes (2008) with a system produced by the Gromacs MD package"), 
						std::string("gmx-morita2008.dat")) { }

		protected:
			//! Selects the waters to be analyzed by loading the entire set of molecules above a particular cutoff location
			void SelectAnalysisWaters ();

	};	// class GMX sfg analyzer

} // namespace morita

#endif
