#include "ambersystem.h"

AmberSystem::AmberSystem (const std::string& prmtop, const std::string& mdcrd, const std::string& mdvel)
// some initialization needs to happen here
: 	_topfile(prmtop),
  _coords(mdcrd, _topfile.NumAtoms()),
  _forces(mdvel, _topfile.NumAtoms())
{
  _atoms = Atom_ptr_vec(_topfile.NumAtoms(), (AtomPtr)NULL);

  // A lot of functionality depends on knowing the system size - so we set it here
  MDSystem::Dimensions (_coords.Dims());

  // Create all the atoms that will be used during analysis
  for (Atom_ptr_vec::iterator it = _atoms.begin(); it != _atoms.end(); it++) {
	*it = new Atom ();
  }

  // and parse all the info out of the topology file into the atoms
  this->_ParseAtomInformation ();
  // and the same for the coordinate and force files
  this->_ParseAtomVectors ();
  // then with all the atoms setup, group them into molecules
  this->_ParseMolecules ();

  return;
}

AmberSystem::~AmberSystem () {
  for (Mol_ptr_vec::iterator it = _mols.begin(); it != _mols.end(); it++) {
	delete *it;
  }
  for (Atom_ptr_vec::iterator it = _atoms.begin(); it != _atoms.end(); it++) {
	delete *it;
  }
}

// While the crdfile holds spatial coordinate information, and the topology file holds atomic information, the data has to be processed into proper atoms in order to play around with them more effectively.
// This function should only be used once when first loading the system up. This info doesn't change during the course of the MD run
void AmberSystem::_ParseAtomInformation () {
  int i = 0;
  for (Atom_ptr_vec::iterator it = _atoms.begin(); it != _atoms.end(); it++) {
	(*it)->Name (_topfile.AtomNames()[i]);
	(*it)->SetAtomProperties();
	(*it)->ID (i);		// set the atom's index number - because we may need to access ordered/list info elsewhere
	++i;
  }

  return;
}

// This function copies the coordinates and force vectors into each of the atoms from the input files
void AmberSystem::_ParseAtomVectors () {

  Atom_it atom_i = _atoms.begin();
  VecR_it coord_i = _coords.begin();
  VecR_it force_i = _forces.begin();
  while (atom_i != _atoms.end()) {
	(*atom_i)->Position (*coord_i);
	if (_forces.Loaded())
	  (*atom_i)->Force (*force_i);

	++atom_i; ++coord_i; ++force_i;
  }

  return;
}

/* the mol pointers in the topology file list the beginning and end atoms of each molecule. This function will group all atoms in a molecule, form a molecule object, and add it to the _mols vector
 */
void AmberSystem::_ParseMolecules () {

  // Here we run through each mol pointer, create a new molecule, and add in the appropriate atoms
  for (int mol = 0; mol < _topfile.NumMols(); mol++) {

	// At each new pointer we create a molecule and start adding in atoms
	// if the molecule is of a specific type for which a class has been created, then let's use that!
	std::string name = _topfile.MolNames()[mol];
	// a topology file doesn't give us the type, so that will be determined when creating the new molecules
	if (name == "no3") {
	  _mols.push_back (new Nitrate());
	}
	else if (name == "hno3") {
	  _mols.push_back (new NitricAcid());
	}
	// we might also have water molecules!
	else if (name == "h2o") {
	  _mols.push_back (new Water());
	}
	else if (name == "so2") {
	  _mols.push_back (new SulfurDioxide());
	}
	else if (name == "dec" || name == "pds") {
	  _mols.push_back (new Decane());
	}
	// otherwise add on a generic molecule
	else {
	  printf ("AmberSystem::_ParseMolecules() -- Couldn't parse the molecule with name '%s'. Don't know what to do with it\n", name.c_str());
	  exit(1);
	}

	_mols[mol]->Name (name);	// set the molecule's residue name

	// Then add all the atoms between the indexes of the molpointers in the topology file
	// MPI had a weird problem **HERE**... fixed June 2007 ~ESS
	int molsize = _topfile.MolSizes()[mol];
	int molpointer = _topfile.MolPointers()[mol];
	for (int atomCount = 0; atomCount < molsize; atomCount++) {
	  // sets the index of the current atom that's being added to the molecule
	  int curAtom = molpointer + atomCount - 1;

	  // Now we're going to bless this new atom with loads of information about itself and its molecule
	  _mols[mol]->MolID (mol);
	  _mols[mol]->AddAtom( _atoms[curAtom] );			// First add the atom into the molecule
	}
  }

  return;
}

void AmberSystem::LoadFirst () {
  _coords.LoadFirst();
  if (_forces.Loaded()) _forces.LoadFirst();
  this->_ParseAtomVectors ();
  return;
}

void AmberSystem::LoadNext () {
  _coords.LoadNext ();							// load up coordinate information from the file
  if (_forces.Loaded()) _forces.LoadNext ();		// also load the force information while we're at it
  this->_ParseAtomVectors ();
  for (Mol_ptr_vec::iterator it = _mols.begin(); it != _mols.end(); it++)
	(*it)->Unset();							// this sets a particular flag on a molecule

  return;
}

void AmberSystem::PrintCRDFile () const {

  for (Atom_it it = _atoms.begin(); it != _atoms.end(); it++) {
	printf ("  % 4.7f  % 4.7f  % 4.7f", (*it)->Position()[x], (*it)->Position()[y], (*it)->Position()[z]);
  }
  printf ("  % 4.7f  % 4.7f  % 4.7f", this->Dims()[x], this->Dims()[y], this->Dims()[z]);

  return;
}
