#include "moritah2o.h"

namespace morita {

  /*	From the morita 2002 water model
		const double BONDLENGTH_EQ	= 0.9575;		// in angstroms
		const double ANGLE_EQ 		= 104.51*M_PI/180.0;		// in radians
		const double CHARGE_H_EQ		= 0.3285;		// charge units (atomic units?)
		const double CHARGE_O_EQ		= -0.6570;
   */

  /* from the equilibrium values of the spcfw model */
/*
  const double BONDLENGTH_EQ	= 1.012 * sfg_units::ANG2BOHR;		// in atomic units
  const double ANGLE_EQ 		= 113.24;		// in degress
  const double CHARGE_H_EQ		= 0.41;		// charge units (atomic units?)
  const double CHARGE_O_EQ		= -0.82;
*/

// using values from the actual simulations
  const double BONDLENGTH_EQ	= 0.99 * sfg_units::ANG2BOHR;		// in atomic units
  const double ANGLE_EQ 		= 105.0;		// in degress
  const double CHARGE_H_EQ		= 0.3285;		// charge units (atomic units?)
  const double CHARGE_O_EQ		= -0.6570;


  // these are the constants from the 2002 Morita work for calculating the dipole moment
  const double C1 =  0.1396,
		C2 = -0.1196,
		C3 = -0.0164,
		C4 = -0.0288,
		C5 = -0.0516,
		C6 =  0.0532,
		C7 = -0.0699,
		C8 =  0.0169,
		C9 =  0.1142;

  const double D1 = 4.6077,
		D2 = 4.8894,
		D3 = 5.5062,
		D4 = 1.6890,
		D5 = 1.6102,
		D6 = 7.3812,
		D7 = 3.4710;


  //
  // Calculates the dipole moment vector of a given water molecule according to the method of Morita/Hynes 2002
  //
  void MoritaH2O::SetDipoleMoment () {

	this->SetAtoms();	// first get the bonds set in the water

	// determine the bondlength displacements
	dR1 = this->OH1()->Magnitude() - BONDLENGTH_EQ;
	dR1 *= sfg_units::ANG2BOHR;
	dR2 = this->OH2()->Magnitude() - BONDLENGTH_EQ;
	dR2 *= sfg_units::ANG2BOHR;

	// angle displacement from the equilibrium
	dA = (acos(this->Angle())*180.0/M_PI) - ANGLE_EQ;

	// symmetry-adapted coordinates
	X1 = dR1 + dR2;
	X2 = dR1 - dR2;

	// determine the charge from the fitting parameters
	qO = CHARGE_O_EQ + C1*X1 + C2*dA + C3*X1*X1 + C4*X2*X2 + C5*dA*dA + C6*X1*dA;
	double dq = C7*X2 + C8*X1*X2 + C9*X2*dA; // originally delta q1 - delta q2
	qH1 = (dq-qO)/2.0;
	qH2 = (-dq-qO)/2.0;

	// the dipole vector is thus classically determined (summing position * charge)
	this->_dipole.Zero();
	this->_dipole += (_o->Position() * qO);
	this->_dipole += (_h1->Position() * qH1);
	this->_dipole += (_h2->Position() * qH2);

  } // Set dipole moment




  void MoritaH2O::SetPolarizability () {
	_alpha1.Zero();

	_alpha1(0,0) = D1+D4*dR1;
	_alpha1(1,1) = D2+D5*dR1;
	_alpha1(2,2) = D3+D6*dR1+D7*dR1*dR1;

	_alpha2(0,0) = D1+D4*dR2;
	_alpha2(1,1) = D2+D5*dR2;
	_alpha2(2,2) = D3+D6*dR2+D7*dR2*dR2;

	_alpha.setZero();

	// These rotations will produce a polarizability tensor that is formed 
	// in the space-fixed frame (instead of the local or molecular frames
	// as written in the paper).
	this->DCMToLabMorita(1);
	_alpha = _alpha + (this->_DCM.transpose() * _alpha1 * this->_DCM);

	this->DCMToLabMorita(2);
	_alpha = _alpha + (this->_DCM.transpose() * _alpha2 * this->_DCM);

  }	// Set Polarizability
  


  DipoleFieldTensor::DipoleFieldTensor (const MoritaH2O_ptr wat1, const MoritaH2O_ptr wat2)
  {
	VecR r = MDSystem::Distance(wat1->GetAtom(Atom::O), wat2->GetAtom(Atom::O));
	// work in atomic units (au)
	double distance = r.Magnitude() * sfg_units::ANG2BOHR;
	double ir3 = 1.0/pow(distance,3.0);
	double ir5 = 3.0/pow(distance,5.0);

	// calculate T as in eq. 10 of the Morita/Hynes 2002 paper
	for (unsigned i = 0; i < 3; i++) {
	  for (unsigned j = 0; j < 3; j++) {
		this->operator()(i,j) = ir3 - ir5*r[i]*r[j];
	  }
	}
  } //dipole field tensor c-tor

} // namespace morita
