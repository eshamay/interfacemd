#include "h2o.h"

int Water::numWaters = 0;

#ifdef H2O_DIPOLE_PARM
WaterDipoleParms Water::_dipparms ("dipoleparm.dat");
#endif

Water::Water ()
{
  this->Rename("h2o");
  ++numWaters;
}

Water::~Water () {
  numWaters--;
}

Water::Water (const Molecule& mol) : Molecule(mol) {
  this->Rename("h2o");
  ++numWaters;
}

Water::Water (const MolPtr& mol) : Molecule(*mol) {
  this->Rename("h2o");
  ++numWaters;
}

void Water::SetBondLengths () {
  this->_oh1 = this->_h1->Position() - this->_o->Position();
  this->_oh2 = this->_h2->Position() - this->_o->Position();
  return;
}

void Water::SetAtoms () {

  // first let's grab pointers to the three atoms and give them reasonable names
  this->_h1 = (Atom *)NULL; this->_h2 = (Atom *)NULL;

  for (Atom_it it = _atoms.begin(); it != _atoms.end(); it++) {
	if ((*it)->Name().find("O") != std::string::npos)
	  this->_o = *it;
	if ((*it)->Name().find("H") != std::string::npos) {
	  if (_h1 == (Atom *)NULL)
		this->_h1 = *it;
	  else
		this->_h2 = *it;
	}

	if (*it == (Atom *)NULL) {
	  std::cout << "problem setting the water atoms! Water::SetAtoms()" << std::endl;
	}
  }

  // we can calculate the two O-H vectors
  this->SetBondLengths ();

  _set = true;
  return;
}

// flip the water about a given plane (perpendicular to the given axis) running through the oxygen
void Water::Flip (const coord axis) {

  this->SetAtoms();

  double center = _o->Position()[axis];
  double distance1 = _oh1[axis];
  double distance2 = _oh2[axis];

  VecR offset1, offset2;

  offset1.Set(axis, center - distance1);
  offset2.Set(axis, center - distance2);

  _h1->Position(axis, center - distance1);
  _h2->Position(axis, center - distance2);

  return;
}

VecR Water::Bisector () {

  this->SetAtoms();

  VecR bisector = _oh1.Unit() + _oh2.Unit();

  return bisector.Unit();
}

/* The molecular axes are defined as follows:
 * 	The Z-axis lies along one of the OH bonds (specified by the index passed)
 * 	The Y-axis is normal to the plane of the water molecule
 * 	The X-axis is then normal to the Z and Y axes, where the 2nd OH bond lies in the positive x-direction
 */
void Water::SetMoritaAxes (const int Zbond) {

  this->SetAtoms();

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


  // the X-axis is just the cross product of the other two
  _x = (_y % _z).Unit();

  return;
}

/* Another type of molecular axes (non-morita, useful for molecular symmetry/order parameter stuff) is as follows:
 * z-axis = the negative of the C2V (bisector) axis
 * y-axis = perpendicular to the plane of the molecule
 * x-axis = y % z
 */
void Water::SetOrderAxes () {

  this->SetAtoms ();

  // the z-axis is the negative of the C2V axis - so find the bisector and set the vector pointing towards the O (just like the dipole vector)
  _z = this->Bisector() * (-1.0);

  // the y-axis points perpendicular to the plane of the molecule. This can be found from the cross product of the two OH vectors
  _y = (_oh1 % _oh2).Unit();

  // and the x-axis is easy
  _x = (_y % _z).Unit();

  return;
}

VecR Water::MolecularAxis () {
  this->SetOrderAxes ();
  return _z;
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
MatR const & Water::DCMToLabMorita (const coord axis, const int bond) {
  this->SetMoritaAxes (bond);

  this->DCMToLab (axis);

  return _DCM;
}

// the alternative is to use the axes definition based on the molecular z-axis lying on the H2O bisector
MatR const & Water::DCMToLabOrder () {
  this->SetOrderAxes ();

  this->DCMToLab ();

  return _DCM;
}

// this should calculate the Euler Angles to get from the molecular frame to the lab frame
void Water::CalcEulerAngles (const coord axis) {

  // First let's set up the direction cosine matrix. The values of the euler angles come from that.
  // Don't forget to set the molecular axes before using this!
  this->DCMToLab (axis);

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

  EulerMatrix.Set (euler_matrix);

  return;
}
