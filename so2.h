#ifndef SO2_H_
#define SO2_H_

#include "molecule.h"

class SulfurDioxide : public Molecule {

  public:
	SulfurDioxide ();
	SulfurDioxide (const Molecule& mol);
	SulfurDioxide (const MolPtr& mol);

	void SetAtoms ();

  private:
	AtomPtr	_s, _o1, _o2;
	VecR	_so1, _so2;

}; // Sulfur Dioxide


#endif
