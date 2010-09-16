#pragma once
#ifndef CTC_H_
#define CTC_H_

#include "molecule.h"

class CarbonTetrachloride : public Molecule {

  public:
	CarbonTetrachloride ();
	CarbonTetrachloride (const Molecule& mol);
	CarbonTetrachloride (const MolPtr& mol);

	void SetAtoms ();

	//double Angle () const { return (_so1 < _so2); } // cos of the O-S-O angle


  //private:
	//AtomPtr	_s, _o1, _o2;
	//VecR	_so1, _so2;

}; // Sulfur Dioxide


#endif
