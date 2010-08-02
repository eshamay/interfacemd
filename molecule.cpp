#include "molecule.h"

int Molecule::numMolecules = 0;

// A constructor for an empty molecule
Molecule::Molecule () :
  _set(false),
  _mass(0.0),
  _name(""),
  _moltype(Molecule::NO_MOLECULE) {
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
  _moltype(oldMol._moltype),
  _DCM (oldMol._DCM) {
  this->Rename(oldMol.Name());
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
AtomPtr Molecule::operator[] (const std::string& atomname) const {
  Atom_it it = std::find_if( _atoms.begin(), _atoms.end(), md_utility::mem_fun_eq(&Atom::Name, atomname));

  // error checking
  if (it == _atoms.end()) {
	printf ("\nFrom Molecule::operator[]\n\"The atom named '%s' was not found in the following molecule:\"\n", atomname.c_str());
	this->Print();
	exit(1);
  }
  return(*it);
}

// get back the first atom pointer to the atom with the given element
AtomPtr Molecule::operator[] (const Atom::Element_t elmt) const {
  
  Atom_it it = std::find_if(_atoms.begin(), _atoms.end(), md_utility::mem_fun_eq(&Atom::Element, elmt));

  // error checking
  if (it == _atoms.end()) {
	printf ("\nFrom Molecule::operator[]\n\"The atom with the given element type was not found in the following molecule:\"\n");
	this->Print();
	exit(1);
  }
  return(*it);
}

AtomPtr Molecule::GetAtom (const std::string& atomname) const {
  AtomPtr patom = (*this)[atomname];
  return(patom);
}

AtomPtr Molecule::GetAtom (const Atom::Element_t elmt) const {
  AtomPtr patom = (*this)[elmt];
  return(patom);
}

// same as the operator, but without the syntax
void Molecule::AddAtom (AtomPtr const atom) {
  _atoms.push_back (atom);            // add the central-periodic image of the atom
  this->FixAtom (atom);
  return;
}

// Instead of adding an atom to a molecule - fix an atom's info in case it was somehow messed up
void Molecule::FixAtom (AtomPtr atom) {

  // now adjust the center of mass to accomodate for the new atom. Also, adjust the total molecular mass.
  //this->UpdateCenterOfMass ();

  // and some of the atom's properties should be adjusted
  atom->ParentMolecule (this);
  atom->Residue (this->_name);
  atom->MolID (this->_ID);

  return;
}

void Molecule::FixAtoms () {
  for (Atom_it atom = _atoms.begin(); atom != _atoms.end(); atom++) {
	this->FixAtom(*atom);
  }
  this->SetAtoms();
}

void Molecule::AddHydrogen (AtomPtr const atom) {

  this->AddAtom (atom);

  if (atom->Element() != Atom::H) {
	std::cout << "Molecule::AddHydrogen() - Tried adding a hydrogen with a non-hydrogen atom:" << std::endl;
	atom->Print();
	std::cout << "To the following molecule" << std::endl;
	this->Print();
	exit(1);
  }

  // now rename the molecule accordingly
  if (_moltype == Molecule::NO3) { _name = "hno3"; _moltype = Molecule::HNO3; }
  else if (_moltype == Molecule::H2O) { _name = "h3o"; _moltype = Molecule::H3O; }
  else if (_moltype == Molecule::OH) { _name = "h2o"; _moltype = Molecule::H2O; }

  this->FixAtoms();

  return;
}

void Molecule::RemoveAtom (AtomPtr const atom) {

  _atoms.erase(std::find(_atoms.begin(), _atoms.end(), atom));

  //this->UpdateCenterOfMass();

  // note on the atom that it is no longer part of a molecule
  atom->ParentMolecule ((MolPtr) NULL);
  atom->Residue ("");
  atom->MolID (-1);

  // if we happen to be taking off a hydrogen from a molecule... rename it accordingly
  if (atom->Element() == Atom::H) {
	if (_moltype == Molecule::HNO3) { _name = "no3"; _moltype = Molecule::NO3; }
	else if (_moltype == Molecule::H3O) { _name = "h2o"; _moltype = Molecule::H2O; }
	else if (_moltype == Molecule::H2O) { _name = "oh"; _moltype = Molecule::OH; }

	this->FixAtoms();
  }

  return;
}

// rename the molecule, and its atoms.
void Molecule::Rename (const std::string& name) {

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
	_centerofmass += (ri * mi);
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

  axis.normalize();

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

  for (Atom_it it = this->begin(); it != this->end(); it++) {
	(*it)->Shift (shift);
  }

  return;
}

void Molecule::Print () const {

  printf ("Residue = %s  (%d)\tmass = % .3f\t%zu wannier centers\n", _name.c_str(), _ID, _mass, _wanniers.size());
  std::for_each (this->begin(), this->end(), std::mem_fun(&Atom::Print));

  printf ("Wannier) ");
  std::for_each (this->_wanniers.begin(), this->_wanniers.end(), std::mem_fun_ref(&VecR::Print));
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
  _moltype = Molecule::NO_MOLECULE;
}

double Molecule::MinDistance (Molecule& mol) {
  // go through the atoms on each molecule and calculate the distance between them, then return the minimum
  bool first = true;
  double min = 0.0;

  for (Atom_it it = this->begin(); it != this->end(); it++) {
	for (Atom_it jt = mol.begin(); jt != mol.end(); jt++) {

	  double temp = ((*it)->Position() - (*jt)->Position()).Magnitude();
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

// returns the direction cosine matrix to the lab frame from the molecular one
MatR const & Molecule::DCMToLab () {
  // These are the three lab-frame axes
  VecR X = Vector3d::UnitX();
  VecR Y = Vector3d::UnitY();
  VecR Z = Vector3d::UnitZ();

  // Here we'll create the lab-frame rotation matrix to rotate molecular properties into the lab-frame
  _DCM.setZero(3,3);
  _DCM << (_x < X) , (_y < X) , (_z < X),
	   (_x < Y) , (_y < Y) , (_z < Y),
	   (_x < Z) , (_y < Z) , (_z < Z);

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
MolPtr Molecule::Merge (MolPtr mol) {
  printf ("merging two molecules:\n");
  this->Print();
  printf ("mol 2---\n");
  mol->Print();
  // the new molecule's name is yet unknown
  _name = "undefined";
  _moltype = Molecule::NO_MOLECULE; 

  for (Atom_it it = mol->begin(); it != mol->end(); it++) {
	this->AddAtom (*it);
  }

  return (this);
}
