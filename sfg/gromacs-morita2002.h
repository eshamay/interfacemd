#ifndef AMBER_MORITA2002_H_
#define AMBER_MORITA2002_H_

#include "morita2002.h"

namespace morita {

	template <class T>
//	class GMXMorita2008Analysis : public Morita2008LookupAnalysis< gromacs::GMXSystem<T> > {
	class GMXMorita2008Analysis : public Morita2002Analysis< gromacs::GMXSystem<T> > {

		public:
			GMXMorita2008Analysis () :
				Morita2002Analysis< gromacs::GMXSystem<T> >(
						std::string("Analysis of the waters in an MD system using the technique of Morita/Hynes (2008) with a system produced by the Gromacs MD package"), 
						std::string("gmx-morita2008.dat")) { }

		protected:
			//! Selects the waters to be analyzed by loading the entire set of molecules above a particular cutoff location
			void SelectAnalysisWaters ();

			//! sets the dipole moments of all the system waters by means of the wannier localization centers calculated by CP2K
			void SetAnalysisWaterDipoleMoments ();
			//! Calculate the polarizability of the given water molecule
			void SetAnalysisWaterPolarizability ();
	};	// class GMX sfg analyzer

} // namespace morita

#endif
