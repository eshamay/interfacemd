#ifndef MOLECULE_FACTORY_H_
#define MOLECULE_FACTORY_H_

#include "molecule.h"
#include "h2o.h"
#include "oh.h"
#include "h.h"
#include "h3o.h"
#include "hno3.h"
#include "so2.h"
#include "ctc.h"
#include "decane.h"

namespace molecule {

	inline MolPtr MoleculeFactory (const std::string name) {

		MolPtr mol;

		if (name == "h2o" || name == "wat") 
			mol = new Water;
		else if (name == "oh")
			mol = new Hydroxide;
		else if (name == "h3o")
			mol = new Hydronium;
		else if (name == "h")
			mol = new Proton;
		else if (name == "so2")
			mol = new SulfurDioxide;
		//else if (name == "no3")
			//mol = new Nitrate;
		else if (name == "ctc")
			mol = new CarbonTetrachloride;

		else {
			mol = new Molecule;
			//std::cerr << "Couldn't determine the molecule using the given name: " << name << std::endl;
			//exit(1);
		}

		return mol;
	}

}	// namespace molecule

#endif
