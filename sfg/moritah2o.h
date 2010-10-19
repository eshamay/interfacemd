#pragma once
#ifndef MORITAH2O_H_
#define MORITAH2O_H_

#include "../mdsystem.h"
#include "sfgunits.h"

namespace morita {

	// a new water that has all the needed pieces for our calculations
	class MoritaH2O : public Water {
		public:

			MoritaH2O (const Molecule& molecule) 
				: Water(molecule) { return; }
			MoritaH2O (const Molecule * molecule) 
				: Water(*molecule) { return; }

			~MoritaH2O () { }

			void SetMoritaAxes (const int bond = 1);		// Determines the molecular-frame axes (a la Morita&Hynes2000) with one bond on the Z-axis, the other in the positive X direction.

			virtual MatR const & DCMToLab ();
			MatR const & DCMToLabMorita (const int bond = 1);	// get the direction cosine matrix for rotations to the lab frame from the morita-hynes one
			MatR const & DCMToLabOrder ();						// direction cosine matrix using the bisector as the molecular z-axis

			MatR EulerMatrix;					// The euler rotation matrix
			double EulerAngles[3];				// euler angles as defined in "The Raman Effect" Appendix A5 (theta, phi, chi)

			void CalcEulerAngles ();

			void SetBondAngleVars ();
			void SetDipoleMoment ();
			void SetDipoleMoment (const VecR& vec) { _dipole = vec; }
			void CalculateGeometricalPolarizability ();
			void SetPolarizability (const MatR& mat) { _alpha = mat; }

			MatR& Polarizability() { return _alpha; }

		protected:

			MatR _alpha1, _alpha2;	// polarizabilities of the two OH bonds
			MatR _alpha;

			double dR1, dR2, dA;		// displacements of the OH bonds and the HOH angle from their equilibrium values
			double X1, X2;
			double qO, qH1, qH2;		// charges calculated from the MoritaHynes 2002 method
	};	// morita-h2o

	typedef MoritaH2O * MoritaH2O_ptr;
	typedef std::vector<MoritaH2O_ptr> Morita_ptr_vec;
	typedef Morita_ptr_vec::const_iterator Morita_it;



}	// namespace morita

#endif
