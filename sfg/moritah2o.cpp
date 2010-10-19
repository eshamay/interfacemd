#include "moritah2o.h"

namespace morita {

	/*	From the morita 2002 water model */
	/*
		 const double BONDLENGTH_EQ	= 0.9575;		// in angstroms
		 const double ANGLE_EQ 		= 104.51*M_PI/180.0;		// in radians
		 const double CHARGE_H_EQ		= 0.3285;		// charge units (atomic units?)
		 const double CHARGE_O_EQ		= -0.6570;
		 */

	/* from the equilibrium values of the spcfw model */
	const double BONDLENGTH_EQ	= 1.012 * sfg_units::ANG2BOHR;		// in atomic units
	const double ANGLE_EQ 		= 113.24;		// in degress
	const double CHARGE_H_EQ		= 0.41;		// charge units (atomic units?)
	const double CHARGE_O_EQ		= -0.82;

	/*
	// using values from the actual cp2k simulations
	const double BONDLENGTH_EQ	= 0.99 * sfg_units::ANG2BOHR;		// in atomic units
	const double ANGLE_EQ 		= 105.0;		// in degress
	const double CHARGE_H_EQ		= 0.3285;		// charge units (atomic units?)
	const double CHARGE_O_EQ		= -0.6570;
	*/


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


	void MoritaH2O::SetBondAngleVars () {

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

		return;
	}

	// Calculates the dipole moment vector of a given water molecule according to the method of Morita/Hynes 2002
	void MoritaH2O::SetDipoleMoment () {
		std::cout << "setting dipole moment" << std::endl;
		this->SetBondAngleVars ();

		// determine the charge from the fitting parameters
		qO = CHARGE_O_EQ + C1*X1 + C2*dA + C3*X1*X1 + C4*X2*X2 + C5*dA*dA + C6*X1*dA;
		double dq = C7*X2 + C8*X1*X2 + C9*X2*dA; // originally delta q1 - delta q2
		// from the relationships that qH1+qH2+qO = 0 (neutrality)
		// and qH1-qH2 = dq (from eq 28 of the 2002 paper)
		qH1 = (dq-qO)/2.0;
		qH2 = (-dq-qO)/2.0;

		// the dipole vector is thus classically determined ( sum(position_i * charge_i) )
		this->_dipole.setZero();
		this->_dipole += (_o->Position() * qO);
		this->_dipole += (_h1->Position() * qH1);
		this->_dipole += (_h2->Position() * qH2);

	} // Set dipole moment




	void MoritaH2O::CalculateGeometricalPolarizability () {
		this->SetBondAngleVars ();

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


	/* The molecular axes are defined as follows:
	 * 	The Z-axis lies along one of the OH bonds (specified by the index passed)
	 * 	The Y-axis is normal to the plane of the water molecule
	 * 	The X-axis is then normal to the Z and Y axes, where the 2nd OH bond lies in the positive x-direction
	 */
	void MoritaH2O::SetMoritaAxes (const int Zbond) {

		this->SetAtoms();

		// try one of the oh bonds as the z-axis
		if (Zbond == 1) {
			_z = _oh1.normalized();
			// let's find the y-axis, as it's just the normal to the molecular plane
			_y = (_oh1 % _oh2).normalized();
		}
		else if (Zbond == 2) {
			_z = _oh2.normalized();
			_y = (_oh2 % _oh1).normalized();
		}


		// the X-axis is just the cross product of the other two
		_x = (_y % _z).normalized();

		return;
	}


	// The molecular axes are defined as per Morita&Hynes (2000) where they set one of the OH bonds (oh1) as the molecular z-axis, and the other bond points in the positive x-axis direction. The result is setting DCM as the direction cosine matrix, that, when operating on a vector in the molecular frame will rotate it into lab-frame coordinates
	MatR const & MoritaH2O::DCMToLabMorita (const int bond) {
		this->SetMoritaAxes (bond);

		Molecule::DCMToLab ();

		return _DCM;
	}

	MatR const & MoritaH2O::DCMToLab () {
		this->DCMToLabMorita();
		return _DCM;
	}

	// the alternative is to use the axes definition based on the molecular z-axis lying on the H2O bisector
	MatR const & MoritaH2O::DCMToLabOrder () {
		this->SetOrderAxes ();

		this->DCMToLab ();

		return _DCM;
	}

	// this should calculate the Euler Angles to get from the molecular frame to the lab frame
	void MoritaH2O::CalcEulerAngles () {

		// First let's set up the direction cosine matrix. The values of the euler angles come from that.
		// Don't forget to set the molecular axes before using this!
		this->DCMToLab ();

		// here is the direct calculation of the euler angles from the direction cosine matrix. This method comes from wikipedia of all places :)
		double x3 = _DCM(0,2);
		double y3 = _DCM(1,2);
		double z1 = _DCM(2,0);
		double z2 = _DCM(2,1);
		double z3 = _DCM(2,2);

		/* If all three axes in the molecular (xyz) and lab (XYZ) frames are aligned, then the euler rotations work by rotating about the body-fixed
		 * axes as follows based on the ZXZ convention:
		 * First a rotation of alpha about the z-axis.
		 * Second a rotation of beta about the x-axis. This is also known as the "tilt" angle because it is the angle between the z and Z axes.
		 * Lastly a rotation of gamma about the body-fixed z-axis. This is the "twist" angle.
		 */
		double alpha = atan2(x3,-y3);
		double beta = acos(z3);
		//double beta = atan2(sqrt(z1*z1+z2*z2), z3);
		double gamma = atan2(z1,z2);

		// alpha is ranged in [-pi/2, pi/2]
		//alpha = fmod(alpha, M_PI/2.0);

		// beta is ranged in [0, pi]
		//beta = fmod(beta, M_PI);

		// because of water's biaxial symmetry, gamma ranges [-pi/2,pi/2]
		//gamma = fmod(gamma, M_PI/2.0);

		//printf ("% 10.4f% 10.4f% 10.4f\n", theta, phi, chi);

		//alpha = -alpha;
		EulerAngles[0] = alpha;
		EulerAngles[1] = beta;
		EulerAngles[2] = gamma;

		double sa = sin(alpha);
		double sb = sin(beta);
		double sg = sin(gamma);
		double ca = cos(alpha);
		double cb = cos(beta);
		double cg = cos(gamma);

		// and now set up the euler rotation matrix according to the ZXZ convention for euler rotations
		double euler_matrix[9] = {
			ca*cg-sa*cb*sg, 	sa*cg+ca*cb*sg, 	sb*sg,
			-ca*sg-sa*cb*cg, 	-sa*sg+ca*cb*cg, 	sb*cg,
			sb*sa, 				-sb*ca,				cb
		};

		MatR t (euler_matrix);
		EulerMatrix = t.transpose();

		return;
	}


} // namespace morita
