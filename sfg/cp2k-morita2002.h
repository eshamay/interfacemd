#pragma once
#ifndef CP2K_MORITA2002_H_
#define CP2K_MORITA2002_H_

#include "morita2002.h"

/*!
 * A derived Morita2002 analysis routine. This performs the SFG spectral analysis (J. Phys. Chem. B 2002, 106, 673-685) on an MD data set produced by the CP2K (Quickstep) ab initio (DFT) molecular dynamics suite (http://cp2k.berlios.de/). This analysis assumes that the input MD data has been calculated along with the corresponding wannier localization centers in order to reproduce molecular and system dipole moments.
 */
namespace morita {

  class CP2KMorita2002Analysis : public Morita2002Analysis<XYZSystem> {
	public:

	  CP2KMorita2002Analysis (std::string filename) : 
		Morita2002Analysis<XYZSystem>(
			std::string("Analysis of the waters in an MD system using the technique of Morita/Hynes (2002) with a system produced by the CP2K package"),
			filename) { }

	protected:
	  //! Selects the waters to be analyzed by loading the entire set - doesn't discriminate on location because there are so few water molecules
	  void SelectAnalysisWaters ();
	  //! sets the dipole moments of all the system waters by means of the wannier localization centers calculated by CP2K
	  void SetAnalysisWaterDipoleMoments ();

  };	// class CP2K sfg analyzer

} // namespace morita

#endif

