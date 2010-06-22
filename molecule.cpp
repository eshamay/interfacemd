#include "molecule.h"

int Molecule::numMolecules = 0;

// A constructor for an empty molecule
Molecule::Molecule () :
	_set(false),
	_mass(0.0),
	_name("")
{
	++numMolecules;
}

// a copy constructor to do a deep copy of a molecule instead of just referencing a pre-existing one.
Molecule::Molecule (const Molecule& oldMol) :
	_atoms(oldMol._atoms),
	_wanniers(oldMol._wanniers),
	_set(false),
	_centerofmass(oldMol._centerofmass),
	_mass(oldMol._mass),
	_name(oldMol._name),
	_ID (oldMol._ID),
	_DCM (oldMol._DCM)
{
	++numMolecules;
}

Molecule::~Molecule () {
	--numMolecules;
}

// Invert the molecule through a point in space. The point is specified by a VecR.
/* The trick here is to shift the entire system such that the point of inversion is at the origin. The invert each coordinate, and return the system back to it's original location.
*/
/*
void Molecule::Invert (VecR origin) {

	for (int i=0; i < _atoms.size(); i++) {
		_atoms[i].Shift (origin * -1.0);
		_atoms[i].Position ( _atoms[i].Position() * -1.0);
		_atoms[i].Shift (origin);
	}
}
*/


/*
	NOTE::: This is broken - the center of mass is only calculated when the molecule is first formed, and not updated at each timestep - fix this!
// For adding another atom into the molecule
int Molecule::operator+= (Atom *atom) {
	this->AddAtom (atom);
return (_atoms.size());
}
*/

// get back the atom pointer to the atom with the given name
Atom * Molecule::operator[] (const std::string atomname) const {

  Atom_it it = std::find_if(_atoms.begin(), _atoms.end(), std::bind2nd(Atom::NameIs_p(), atomname));

  // error checking
  if (it == _atoms.end()) {
	printf ("\nFrom Molecule::operator[]\n\"The atom named '%s' was not found in the following molecule:\"\n", atomname.c_str());
	this->Print();
	exit(1);
  }
	// note here that the atom may not even exist in the molecule! If you try to operate on that pointer something awefully bad could happen. Before using this method for retrieving atoms, make sure that the calling function checks that it knows what to do with a NULL pointer
	
	return(*it);
}

Atom * Molecule::GetAtom (const string atomname) const {
  Atom * patom = (*this)[atomname];
  return(patom);
}

// same as the operator, but without the syntax
void Molecule::AddAtom (Atom * const atom) {
  _atoms.push_back (atom);            // add the central-periodic image of the atom
  this->FixAtom (atom);
  return;
}

// Instead of adding an atom to a molecule - fix an atom's info in case it was somehow messed up
void Molecule::FixAtom (Atom * const atom) {

  // now adjust the center of mass to accomodate for the new atom. Also, adjust the total molecular mass.
  this->UpdateCenterOfMass ();

  // and some of the atom's properties should be adjusted
  atom->ParentMolecule (this);
  atom->Residue (this->_name);
  atom->MolID (this->_ID);

  return;
}

void Molecule::AddHydrogen (Atom * const atom) {

  this->AddAtom (atom);

  if (atom->Name().find("H") == string::npos) {
	std::cout << "Molecule::AddHydrogen() - Tried adding a hydrogen with a non-hydrogen atom:" << std::endl;
	atom->Print();
	exit(1);
  }

  // now rename the molecule accordingly
  if (_name == "no3") this->Rename ("hno3");
  else if (_name == "h2o") this->Rename ("h3o");
  else if (_name == "oh") this->Rename ("h2o");

  return;
}

void Molecule::RemoveAtom (Atom * const atom) {

  _atoms.erase(std::find(_atoms.begin(), _atoms.end(), atom));

  this->UpdateCenterOfMass();

  // note on the atom that it is no longer part of a molecule
  atom->ParentMolecule ((Molecule *)NULL);
  atom->Residue ("");
  atom->MolID (-1);

  // if we happen to be taking off a hydrogen from a molecule... rename it accordingly
  if (atom->Name() == "H") {
	if (_name == "hno3") this->Rename ("no3");
	else if (_name == "h3o") this->Rename ("h2o");
	else if (_name == "h2o") this->Rename ("oh");
  }

  return;
}

// rename the molecule, and its atoms.
void Molecule::Rename (const std::string name) {

  this->_name = name;
  for (Atom_it atom = _atoms.begin(); atom != _atoms.end(); atom++) {
	(*atom)->Residue(name);
  }

  return;
}

// recalculate the center of mass if coordinates are updated
VecR Molecule::UpdateCenterOfMass () {
  // first zero it out
  _centerofmass.Zero();
  _mass = 0.0;

  // then run through each atom's coords and add in its contribution.
  for (Atom_it atom = _atoms.begin(); atom != _atoms.end(); atom++) {
	VecR ri = (*atom)->Position();
	double mi = (*atom)->Mass();
	_mass += mi;
	_centerofmass += (ri) * mi;
  }
  _centerofmass /= _mass;

  return(_centerofmass);
}


/*
// join two molecules into one
int Molecule::operator+= (Molecule& mol) {
for (int i = 0; i < mol.size(); i++) {
(*this) += (mol[i]);
}

return (_atoms.size());
}
*/

/* This will reflect (invert) about a given axis */

void Molecule::Reflect (coord const axis, double const plane) {

  for (Atom_it atom = _atoms.begin(); atom != _atoms.end(); atom++) {
	VecR pos ((*atom)->Position());
	VecR force ((*atom)->Force());
	(*atom)->Position (axis, 2.0*plane - pos[axis]);
	(*atom)->Force (axis, -force[axis]);
  }
}

/*
// pbase_Plane(): get base of perpendicular from point to a plane
//    Input:  P = a 3D point
//            PL = a plane with point V0 and normal n
//    Output: *B = base point on PL of perpendicular from P
//    Return: the distance from P to the plane PL
pbase_Plane( Point point, Plane plane, Point* B)
float    sb, sn, sd;

sn = -dot( plane.n, (point - plane.V0));
sd = dot(plane.n, plane.n);
sb = sn / sd;

 *B = point + sb * plane.n;
 return d(point, *B);

*/


void Molecule::Rotate (VecR& origin, VecR& axis, double angle) {

  double nx = 0.0, ny = 0.0, nz = 0.0;

  axis = axis.Unit();


  for (Atom_it atom = _atoms.begin(); atom != _atoms.end(); atom++) {
	// first move all the atoms to the origin of the rotation axis
	(*atom)->Shift (origin * -1.0);

	// then apply the rotation as per wikipedia's 'Rotation Matrix' entry
	nx = (cos(angle)+(1-cos(angle))*pow(axis[x],2))*(*atom)->X() +
	  ((1-cos(angle))*axis[x]*axis[y]-sin(angle)*axis[z])*(*atom)->Y() +
	  ((1-cos(angle))*axis[x]*axis[z]+sin(angle)*axis[y])*(*atom)->Z();

	ny = ((1-cos(angle))*axis[y]*axis[x]+sin(angle)*axis[z])*(*atom)->X() +
	  (cos(angle)+(1-cos(angle))*pow(axis[y],2))*(*atom)->Y() +
	  ((1-cos(angle))*axis[y]*axis[z]-sin(angle)*axis[x])*(*atom)->Z();

	nz = ((1-cos(angle))*axis[z]*axis[x]-sin(angle)*axis[y])*(*atom)->X() +
	  ((1-cos(angle))*axis[z]*axis[y]+sin(angle)*axis[x])*(*atom)->Y() +
	  (cos(angle)+(1-cos(angle))*pow(axis[z],2))*(*atom)->Z();

	// having calc'd the new positions about the axis, move the atom
	(*atom)->Position (VecR (nx, ny, nz));

	// and then relocate it back to the original origin
	(*atom)->Shift (origin);
  }
}

void Molecule::Shift (VecR& shift) {

  RUN (_atoms) {
	_atoms[i]->Shift (shift);
  }

  return;
}

/*
// note: this is going to CREATE AN ARRAY - make sure to DEALLOCATE (i.e. delete) it after use!!
double * Molecule::DPositions () const {
double * positions;

// Each atom has a position vector consisting of three components. So let's allocate space for all these vectors
positions = new double[_atoms.size() * 3];

// now we assign the value from the atomic positions
RUN (_atoms) {
double * pos = _atoms[i]->DPosition();	// this is the position vector of the atom

for (int j = 0; j < 3; j++) {
positions[i * 3 + j] = pos[j];	// and now we set each element in the output array
}
}

return (positions);
}
*/

/*
// Just like the DPositions (), make sure to deallocate the array we create here!
double * Molecule::DForces () const {
double * forces;

// Each atom has a position vector consisting of three components. So let's allocate space for all these vectors
forces = new double[_atoms.size() * 3];

// now we assign the value from the atomic positions
RUN (_atoms) {
double * force = _atoms[i]->DForce();	// this is the position vector of the atom

for (int j = 0; j < 3; j++) {
forces[i * 3 + j] = force[j];	// and now we set each element in the output array
}
}

return (forces);
}
*/

void Molecule::Print () const {

  printf ("Residue = %s  (%d)\tmass = % .3f\n", _name.c_str(), _ID, _mass);

  RUN (_atoms) {
	_atoms[i]->Print();
  }
  //printf ("wan)\t%d\n", _wanniers.size());

  return;
}

void Molecule::clear () {
  _atoms.clear();
  _wanniers.clear();
  _mass = 0.0;
  _centerofmass.Zero();
  //_dipole.Zero();
  _name = "";
  _ID = 0;
}

// This should calculate the dipole of a molecule given that we've already generated the wannier centers
VecR Molecule::CalcDipole () {

  this->UpdateCenterOfMass();

  _dipole.Zero();

  // the dipole is just a sum of the position vectors multiplied by the charges (classical treatment)
  for (Atom_it it = _atoms.begin(); it != _atoms.end(); it++) {
	VecR r = (*it)->Position() - _centerofmass;
	_dipole += r * (*it)->Charge();
  }

  // wannier centers have a charge of -2
  RUN (_wanniers) {
	_dipole -= (_wanniers[i] - _centerofmass) * 2.0;
  }

  return (_dipole);
}

double Molecule::MinDistance (Molecule& mol) {

  // go through the atoms on each molecule and calculate the distance between them, then return the minimum
  bool first = true;
  double min = 0.0;

  RUN (_atoms) {
	RUN2 (mol.Atoms()) {

	  Atom * atom1 = _atoms[i];
	  Atom * atom2 = mol.Atoms(j);

	  double temp = (atom1->Position() - atom2->Position()).Magnitude();
	  //printf ("%f\n", temp);
	  if (first) {
		first = false;
		min = temp;
	  }

	  if (temp < min) {
		min = temp;
	  }
	}
  }

  return (min);
}

// This builds a rotation matrix to rotate from the Lab-frame to the body-fixed frame coordinates
void Molecule::RotateToMol (double vector[3]) const {

  // here's the lab-frame coordinates
  VecR X (1.0, 0.0, 0.0);
  VecR Y (0.0, 1.0, 0.0);
  VecR Z (0.0, 0.0, 1.0);

  // now we build our direction cosine matrix (each element is the cosine of the angle between two axes
  // The molecular-frame axes are already known as _x, _y, and _z
  double rotation[3][3] = {
	{ _x < X, _x < Y, _x < Z},
	{ _y < X, _y < Y, _y < Z},
	{ _z < X, _z < Y, _z < Z}
  };

  this->RotateVector(rotation, vector);

  return;
}

void Molecule::RotateToMol (double matrix[][3]) const {

  // here's the lab-frame coordinates
  VecR X (1.0, 0.0, 0.0);
  VecR Y (0.0, 1.0, 0.0);
  VecR Z (0.0, 0.0, 1.0);

  // now we build our direction cosine matrix (each element is the cosine of the angle between two axes
  // The molecular-frame axes are already known as _x, _y, and _z
  double rotation[3][3] = {
	{ _x < X, _x < Y, _x < Z},
	{ _y < X, _y < Y, _y < Z},
	{ _z < X, _z < Y, _z < Z}
  };

  this->RotateMatrix(rotation, matrix);

  return;
}

// This builds a rotation matrix to rotate from the Lab-frame to the body-fixed frame coordinates
void Molecule::RotateToLab (double vector[3]) const {

  // here's the lab-frame coordinates
  VecR X (1.0, 0.0, 0.0);
  VecR Y (0.0, 1.0, 0.0);
  VecR Z (0.0, 0.0, 1.0);

  // now we build our direction cosine matrix [transpose of the reverse operation - RotateToMol]
  // (each element is the cosine of the angle between two axes)
  // The molecular-frame axes are already known as _x, _y, and _z
  double rotation[3][3] = {
	{ _x < X, _y < X, _z < X},
	{ _x < Y, _y < Y, _z < Y},
	{ _x < Z, _y < Z, _z < Z}
  };

  this->RotateVector(rotation, vector);

  return;
}

void Molecule::RotateToLab (double matrix[][3]) const {

  // here's the lab-frame coordinates
  VecR X (1.0, 0.0, 0.0);
  VecR Y (0.0, 1.0, 0.0);
  VecR Z (0.0, 0.0, 1.0);

  // now we build our direction cosine matrix [transpose of the reverse operation - RotateToMol]
  // (each element is the cosine of the angle between two axes)
  // The molecular-frame axes are already known as _x, _y, and _z
  double rotation[3][3] = {
	{ _x < X, _y < X, _z < X},
	{ _x < Y, _y < Y, _z < Y},
	{ _x < Z, _y < Z, _z < Z}
  };

  this->RotateMatrix(rotation, matrix);

  return;
}

void Molecule::RotateVector (double rotation[][3], double vector[3]) const {

  double temp[3];		// a temp holder for the rotated vector

  for (int i = 0; i < 3; i++) {
	temp[i] = rotation[i][0]*vector[0] + rotation[i][1]*vector[1] + rotation[i][2]*vector[2];
  }

  // and lastly apply the changes to the input vector
  for (int i = 0; i < 3; i++) {
	vector[i] = temp[i];
  }

  return;
}

void Molecule::RotateMatrix (double rotation[][3], double matrix[][3]) const {

  // apply the matrix rotation to the given matrix (probably something like the polarizability derivative matrix)

  double temp[3][3];	// this is a temp output matrix to work with

  for (int i = 0; i < 3; i++) {
	for (int j = 0; j < 3; j++) {
	  temp[i][j] = rotation[i][0]*matrix[0][j] + rotation[i][1]*matrix[1][j] + rotation[i][2]*matrix[2][j];
	}
  }

  // and lastly apply the changes to the input matrix
  for (int i = 0; i < 3; i++) {
	for (int j = 0; j < 3; j++) {
	  matrix[i][j] = temp[i][j];
	}
  }

  return;
}


// returns the direction cosine matrix to the lab frame from the molecular one
MatR const & Molecule::DCMToLab (const coord axis) {
  // These are the three lab-frame axes
  VecR X, Y, Z;

  if (axis == z) {
	X.Set (1.,0.,0.);
	Y.Set (0.,1.,0.);
	Z.Set (0.,0.,1.);
  }

  // but if the system is funky and we want to use different lab-frame coords, perhaps treating the Y-axis as the primary axis, then this will work
  // This is just changing the 'primary' axis used to define the 'tilt' angle when calculating Euler angles
  if (axis == y) {
	X.Set (0.,1.,0.);
	Y.Set (0.,0.,1.);
	Z.Set (1.,0.,0.);
  }

  // Here we'll create the lab-frame rotation matrix to rotate molecular properties into the lab-frame
  double rotation_data[9] = {	_x<X, _y<X, _z<X,
	_x<Y, _y<Y, _z<Y,
	_x<Z, _y<Z, _z<Z   };
  _DCM.Set(rotation_data);

  return _DCM;
}

/* The technique for axes rotations is described well in Zare's book on angular momentum in the appendix of direction cosines.
 * _eulerangles is set up as [theta, phi, chi]
 */
void Molecule::_FindEulerAngles () {

  // Theta - The angle between the reference and derived z-axes
  _eulerangles[0] = acos(_z[z]);

  // Phi, is then found similarly with x and y components from other axes
  _eulerangles[1] = atan2(_y[z], _x[z]);

  /*
  // Chi is not so obvious. We have to actually rotate the axes to find that last angle
  double rotate[3][3];	// the rotation matrix
  double angles[3] = {_eulerangles[0], _eulerangles[1], 0.0};
  VecR rotatedZ = _oh2; // the rotated molecular Z-axis

  this->RotateToLabLabFrame (rotate, angles);
  this->RotateVector (rotate, rotatedZ.Coords());

  // Chi is derived from a combination of the x and y components of the rotated z-axis in the reference frame
  _eulerangles[2] = atan2(rotatedZ[y], rotatedZ[x]);	// Chi/Psi (depending on your choice of greek)
  */

  // Chi is derived from a combination of the x and y components of the rotated z-axis in the reference frame
  _eulerangles[2] = atan2(-_z[y], _z[x]);

  return;
}

/*
   void Molecule::ClearHBonds () {
   RUN (_atoms) {
   _atoms[i]->ClearHBonds();
   }

   return;
   }

// generate a list of all the atoms that are hydrogen-bonding to atoms within this molecule
std::vector<Atom *> Molecule::HBonds () const {

std::vector<Atom *> atoms, atomBonds;
atoms.clear();

RUN (_atoms) {
atomBonds = _atoms[i]->HBonds();

RUN2 (atomBonds) {
atoms.push_back (atomBonds[j]);
}
}

return (atoms);
}
*/

// if given a 2nd molecule, this will merge the current and the new molecules into one larger molecule.
Molecule * Molecule::Merge (Molecule * mol) {

  printf ("merging two molecules:\n");
  this->Print();
  printf ("mol 2---\n");
  mol->Print();
  // the new molecule's name is yet unknown
  _name = "undefined";

  RUN (mol->Atoms()) {
	this->AddAtom (mol->Atoms(i));
  }

  return (this);
}
