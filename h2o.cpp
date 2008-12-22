#include "h2o.h"

int Water::numWaters = 0;

#ifdef H2O_DIPOLE_PARM
WaterDipoleParms Water::_dipparms ("dipoleparm.dat");
#endif

Water::Water () {
	_centerofmass = VecR ();
	_mass = 0.0;
	_name = "h2o";
	_set = false;
	++numWaters;
}

Water::~Water () {
	numWaters--;
}

Water::Water (const Molecule& molecule) {
	// copy over the old information
	_centerofmass = molecule.CenterOfMass ();
	_mass = molecule.Mass ();
	_name = molecule.Name ();

	// now run through and make copies of all the atoms
	RUN (molecule.Atoms()) {
		// this calls up the copy constructor and forms a real copy of the template atom
		Atom *newatom = new Atom(*molecule[i]);
		_atoms.push_back (newatom);
	}
	_set = false;
	++numWaters;
}
	
void Water::FindOHBonds () {
	// first let's grab pointers to the three atoms and give them reasonable names
	_h1 = (Atom *)NULL; _h2 = (Atom *)NULL;

	RUN (_atoms) {
		if (_atoms[i]->Name().find("O") != string::npos)
			_o = _atoms[i];

		if (_atoms[i]->Name().find("H") != string::npos) {
			if (_h1 == (Atom *)NULL)
				_h1 = _atoms[i];
			else
				_h2 = _atoms[i];
		}
	}
	
	// we can calculate the two O-H vectors
	_oh1 = _o->Position().MinVector(_h1->Position(), Atom::Size());
	_oh2 = _o->Position().MinVector(_h2->Position(), Atom::Size());

	_set = true;

return;
}

VecR Water::Bisector () {

	if (!_set) this->FindOHBonds();

	VecR bisector = _oh1.Unit() + _oh2.Unit();

return bisector.Unit();
}

void Water::CalcDipole () {
	
	// this plan goes along the same method of Morita/Hynes (J. Phys. Chem. B, Vol. 106, No. 3, 2002) for calculating the dipole moment. It's based on a parameterization of the water partial charges
	if (!_set) this->FindOHBonds();

	// first we calculate the displacements of the OH bond lengths and the water angle
	double dR1 = _oh1.Magnitude() - R_eq;
	double dR2 = _oh2.Magnitude() - R_eq;
	double dTheta = acos(_oh1 < _oh2) * 180.0 / M_PI - Theta_eq;

	//printf ("%f\t%f\t%f\n", dR1, dR2, dTheta);
	//printf ("%f\t%f\t%f\n", _oh1.Magnitude(), _oh2.Magnitude(), acos(_oh1 < _oh2) * 180.0/M_PI);

	// then calculate the symmetry-adapted products
	double X1 = dR1 + dR2;
	double X2 = dR1 - dR2;
	double X3 = dTheta;

	// these are the constants from the work
	double C1 = 0.1396, C2 = -0.1196, C3 = -0.0164, C4 = -0.0288, C5 = -0.0516, C6 = 0.0532, C7 = -0.0699, C8 = 0.0169, C9 = 0.1142;

	double alpha = C7*X2 + C8*X1*X2 + C9*X2*X3;		// originally delta q1 - delta q2
	double dqO = C1*X1 + C2*X3 + C3*X1*X1 + C4*X2*X2 + C5*X3*X3 + C6*X1*X3;	// originally delta q0
	
	// and calculating the partial charges...
	double qO = qO_eq + dqO;		// partial charge on the oxygen
	double q1 = (alpha - qO) / 2.0; 
	double q2 = (-alpha - qO) / 2.0;

	//printf ("%f\t%f\t%f\t%f\t%f\n", alpha, beta, q0, q1, q2);
	// now the dipole moment is calculated classically
	_dipole = _o->Position()*qO + _h1->Position()*q1 + _h2->Position()*q2;
	
return; 
}

/* The molecular axes are defined as follows:
 * 	The Z-axis lies along one of the OH bonds (specified by the index passed)
 * 	The Y-axis is normal to the plane of the water molecule
 * 	The X-axis is then normal to the Z and Y axes, where the 2nd OH bond lies in the positive x-direction
 */
void Water::SetMoritaAxes (const int Zbond) {
	
	if (!_set) this->FindOHBonds();

	// try one of the oh bonds as the z-axis
	if (Zbond == 1) {
		_z = _oh1.Unit();
		// let's find the y-axis, as it's just the normal to the molecular plane
		_y = (_oh1 % _oh2).Unit();
	}
	else if (Zbond == 2) {
		_z = _oh2.Unit();
		_y = (_oh2 % _oh1).Unit();
	}

	// the X-axis is just the cross produce of the other two
	_x = (_y % _z).Unit();

return;
}

/*
#ifdef WATER_POLARIZ
// This calculations of the molecular polarizability is based on the Morita-Hynes 2002 paper where they use ab initio calcs to find fitting parameters to calculate alpha. We calculate a polarizability tensor for both OH bonds, and then rotate each one into the molecular frame, and then sum them to get the molecular polarizability.
void Water::CalcAlpha () {
	
	double 	D1 = 4.6077, 
			D2 = 4.8894,
			D3 = 5.5062,
			D4 = 1.6890,
			D5 = 1.6102,
			D6 = 7.3812,
			D7 = 3.4710;

	this->FindOHBonds();
	
	double oh1 = _oh1.Magnitude();
	double oh2 = _oh2.Magnitude();

	MatR alpha_1; 
	alpha_1.Set (0,0, oh1*D4 + D1);
	alpha_1.Set (1,1, oh1*D5 + D2);
	alpha_1.Set (2,2, oh1*D6 + oh1*oh1*D7 + D3);

	MatR alpha_2; 
	alpha_2.Set (0,0, oh2*D4 + D1);
	alpha_2.Set (1,1, oh2*D5 + D2);
	alpha_2.Set (2,2, oh2*D6 + oh2*oh2*D7 + D3);

	this->SetMoritaAxes (1);
	
	VecR frame[3];
	frame[0] = _x; frame[1] = _y; frame[2] = _z;
	MatR alpha_1_rot = alpha_1.RotateToFrame (frame);

	this->SetMoritaAxes (2);
	frame[0] = _x; frame[1] = _y; frame[2] = _z;
	MatR alpha_2_rot = alpha_2.RotateToFrame (frame);

	_alpha = alpha_1_rot + alpha_2_rot;

return;
}
#endif
*/

// The molecular axes are defined as per Morita&Hynes (2000) where they set one of the OH bonds (oh1) as the molecular z-axis, and the other bond points in the positive x-axis direction. The result is setting DCM as the direction cosine matrix, that, when operating on a vector in the molecular frame will rotate it into lab-frame coordinates
MatR const & Water::DCMToLab (const int bond) {
    this->SetMoritaAxes (bond);

    // We already have the three molecule-frame axes,

    // These are the three lab-frame axes
    VecR X (1,0,0);
    VecR Y (0,1,0);
    VecR Z (0,0,1);

    // Here we'll create the lab-frame rotation matrix to rotate molecular properties into the lab-frame
    double rotation_data[9] = {_x<X, _x<Y, _x<Z, _y<X, _y<Y, _y<Z, _z<X, _z<Y, _z<Z};
    DCM.Set (rotation_data);

	// while we're at it, let's quickly compute the Euler angles

return DCM;
}

// this should calculate the Euler Angles to get from the molecular frame to the lab frame
double * Water::CalcEulerAngles (const int bond) {

	// First let's set up the direction cosine matrix. The values of the euler angles come from that.
	this->DCMToLab (bond);

	// here is the direct calculation of the euler angles from the direction cosine matrix. This method comes from 
	double z1 = DCM.Index(2,0);
	double z2 = DCM.Index(2,1);
	double z3 = DCM.Index(2,2);

	double beta = atan2(sqrt(z1*z1+z2*z2), z3);
	double alpha = atan2(DCM.Index(0,2),-DCM.Index(1,2));
	double gamma = atan2(DCM.Index(2,0),DCM.Index(2,1));

	//printf ("% 10.4f% 10.4f% 10.4f\n", theta, phi, chi);

	//alpha = -alpha;
	EulerAngles[0] = alpha;
	EulerAngles[1] = beta;
	EulerAngles[2] = gamma;

	// and now set up the euler rotation matrix
    double euler_matrix[9] = {
		cos(alpha)*cos(gamma) - sin(alpha)*cos(beta)*sin(gamma),
		sin(alpha)*cos(gamma) + cos(alpha)*cos(beta)*sin(gamma),
		sin(beta)*sin(gamma),
		-cos(alpha)*sin(gamma)-sin(alpha)*cos(beta)*cos(gamma),
		-sin(alpha)*sin(gamma)+cos(alpha)*cos(beta)*cos(gamma),
		sin(beta)*cos(gamma),
		sin(beta)*sin(alpha),
		-sin(beta)*cos(alpha),
		cos(beta)
	};

    EulerMatrix.Set (euler_matrix);

return EulerAngles;
}
	
int Hydroxide::numHydroxides = 0;

Hydroxide::Hydroxide () {
	_centerofmass = VecR ();
	_mass = 0.0;
	_name = "oh";
	_set = false;
	++numHydroxides;
}

Hydroxide::~Hydroxide () {
	--numHydroxides;
}

void Hydroxide::SetAtoms () {
	_o = (*this)["O"];
	_h = (*this)["H"];

	_oh = _h->Position().MinVector(_o->Position(), Atom::Size());
	_set = true;

return;
}

int Hydronium::numHydroniums = 0;

Hydronium::Hydronium () {
	_centerofmass = VecR ();
	_mass = 0.0;
	_name = "h3o";

	_set = false;
	++numHydroniums;
}

Hydronium::~Hydronium () {
	--numHydroniums;
}

void Hydronium::SetAtoms () {
	
		// here's the hydrogen and nitrogen atoms
		_o = (*this)["O"];

		// now set the 3 hydrogens
		std::vector<Atom *> hydrogens;
		RUN (_atoms) {
			if (_atoms[i]->Name() != "H") continue;
			hydrogens.push_back (_atoms[i]);
		}

		_h1 = hydrogens[0];
		_h2 = hydrogens[1];
		_h3 = hydrogens[2];

		// while we're here we may as well also find the N-O bond vectors
		_oh1 = _h1->Position().MinVector(_o->Position(), Atom::Size());
		_oh2 = _h2->Position().MinVector(_o->Position(), Atom::Size());
		_oh3 = _h3->Position().MinVector(_o->Position(), Atom::Size());

		_set = true;

return;
}

//testing dipole
