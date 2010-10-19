#include "moleculefactory.h"

namespace molecule {

	MolPtr MoleculeFactory (const std::string name) {

		MolPtr mol;

		if (name == "h2o" || name == "wat" || name == "SW") 
			mol = new Water;
		else if (name == "oh")
			mol = new Hydroxide;
		else if (name == "h3o")
			mol = new Hydronium;
		else if (name == "h")
			mol = new Proton;
		else if (name == "so2" || name == "sog" || name == "soa")
			mol = new SulfurDioxide;
		//else if (name == "no3")
		//mol = new Nitrate;
		else if (name == "ctc")
			mol = new CarbonTetrachloride;

		else {
			//mol = new Molecule;
			std::cerr << "Couldn't determine the molecule using the given name: " << name << std::endl;
			exit(1);
		}

		return mol;
	}	// MoleculeFactory

}	// namespace molecule
