#include "xyzsystem.h"

XYZSystem::XYZSystem (string filepath, VecR size, string wannierpath) :
  MDSystem (),
  _coords(filepath),
  _wanniers(wannierpath),
  _reparse_limit(1),	// initially set to parse everything everytime
  _reparse_step(0)
{
  MDSystem::Dimensions (size);
  this->LoadFirst();
}



XYZSystem::~XYZSystem () {
  for (Mol_it it = _mols.begin(); it != _mols.end(); it++)
	delete *it;
}



void XYZSystem::_ParseMolecules () {

  /***********************************************************************************
   * This is the top-level parsing routine to give the overall idea of what's going on
   * *********************************************************************************/

  // first things first - we need the interatomic distances!
  graph.UpdateGraph (_atoms);

  // Now let's do some house-cleaning to set us up for working with new molecules
  for (Mol_it it = _mols.begin(); it != _mols.end(); it++) {
	delete *it;				// get rid of all the molecules in memory
  }
  _mols.clear();				// then clear out the molecule list

  // we also have to go through and clear out some info on all the atoms
  // like the parentmolecules and names
  for (Atom_it it = _atoms.begin(); it != _atoms.end(); it++) {
	(*it)->ParentMolecule ( (Molecule *)NULL );
	(*it)->Residue ("");
	(*it)->MolID (-1);
  }

  // track which atoms have already been parsed in the system
  _unparsed.clear();
  std::copy (_atoms.begin(), _atoms.end(), std::back_inserter(_unparsed));

  this->_ParseWaters ();
  this->_ParseNitrates ();
  this->_ParseSulfides ();

  this->_CheckForUnparsedAtoms ();

  return;
}

void XYZSystem::_CheckForUnparsedAtoms () const {
  bool leave = false;
  for (Atom_it it = _unparsed.begin(); it != _unparsed.end(); it++) {
	if (*it != (AtomPtr)NULL) {
	  (*it)->Print();	// show all remaining atoms

	  // and all the atoms to which it is bound (and the distance)
	  std::vector<AtomPtr> bound (graph.BondedAtoms(*it));
	  for (Atom_it jt = bound.begin(); jt != bound.end(); jt++) { 
		cout << "^--~ (" << graph.Distance(*jt, *it) << ")  ";
		(*jt)->Print();
	  }
	  leave = true;
	}
  }
  if (leave) {
	cout << "Found the above atoms unaccounted for in the system! Fix up the parsing routine" << endl;
	exit(1);
  }
  return;
}


/******************************
 * Add in the wannier centers
 *
 * What we'll do here is run through each of the wannier centers and find which molecule they belong to. For each center we'll check its distance against all the atoms and when we find the one it's bound to we'll shove it into the parent molecule.
 * ****************************/
void XYZSystem::_ParseWanniers () {
  // Then we've got to clear out all the wanniers already loaded into the molecules
  for (Mol_it it = _mols.begin(); it != _mols.end(); it++) {
	(*it)->ClearWanniers();
  }

  for (Atom_it O = _atoms.begin(); O != _atoms.end(); O++) {
	// find every oxygen atom
	if ((*O)->Name().find("O") == string::npos && (*O)->Name().find("S") == string::npos) continue;

	MolPtr mol = (*O)->ParentMolecule();

	// and then go through all the wanniers in order to find the ones attached to the oxygen
	for (VecR_it vi = _wanniers.begin(); vi != _wanniers.end(); vi++) {

	  // we find the distance to the oxygen from the wannier center
	  double distance = MDSystem::Distance ((*O)->Position(), *vi).Magnitude();

	  // if it's close enough, then that oxygen gets it
	  if (distance < WANNIER_BOND) {
		mol->AddWannier(*vi);
	  }
	}
  }

return;
}	// Parse Wanniers



void XYZSystem::LoadFirst () {
  _atoms.clear();
  std::copy(_coords.begin(), _coords.end(), std::back_inserter(_atoms));

  this->_ParseMolecules();

  if (_wanniers.Loaded()) {
	this->_ParseWanniers();
  }
}	// Load First 


// watch out - no functionality for wannier centers here (yet)
void XYZSystem::Seek (int step) {
  _coords.Seek (step);
  this->_ParseMolecules();
}	// Seek


void XYZSystem::LoadNext () {
  _coords.LoadNext();

  if (_reparse_step++ == _reparse_limit) {
	this->_ParseMolecules();
	_reparse_step = 0;
  }

  if (_wanniers.Loaded()) {
	_wanniers.LoadNext();
	this->_ParseWanniers();
  }
} // Load Next



// This will calculate the total dipole moment of the system based on atom locations and wannier center positions
// The origin is shifted to the center of the system in order to get closest images (wrapped into the box) of all the atoms/wanniers
VecR XYZSystem::SystemDipole () {

  VecR dipole;
  //VecR center = Atom::Size() * 0.5;	// center of the system

  for (Atom_it it = _atoms.begin(); it != _atoms.end(); it++) {
	//VecR ri = center.MinVector(_atoms[i]->Position(), Atom::Size());
	VecR ri = (*it)->Position();
	dipole += ri * (*it)->Charge();
  }

  for (VecR_it it = _wanniers.begin(); it != _wanniers.end(); it++) {
	//VecR ri = center.MinVector(_wanniers[i], Atom::Size());
	dipole -= (*it) * 2.0;
  }

  return (dipole);
}

void XYZSystem::_ParseWaters () {

  /*******************
   * Processing Waters
   *******************/

  // let's go through the system and find all the waters and form molecules out of them
  for (Atom_it it = _atoms.begin(); it != _atoms.end(); it++) {

	// every water has an oxygen atom, right? So let's find those
	AtomPtr O = *it;
	if (O->Name().find("O") == string::npos) continue;

	// The bondgraph provides us with all the covalently bound H's
	Atom_ptr_vec atoms = graph.BondedAtoms (O, covalent);
	Atom_ptr_vec Hs = graph.BondedAtoms (O, covalent, "H");


	// if we pick up an OH group that is part of a larger molecule (i.e. nitric acid) then it will be processed as a hydroxide...
	// So if the only type of atom attached is a hydrogen, we have some form of water (OH, H2O, H3O)
	if (Hs.size() != atoms.size()) continue;

	int molIndex = (int)_mols.size();	// set the molecule index

	int num_H = (int)atoms.size();
	// we may be dealing with a hydroxide ion
	if (num_H == 1) {
	  _mols.push_back (new Hydroxide ());
	}

	// otherwise we have a full molecule
	else if (num_H == 2) {
	  _mols.push_back (new Water ());
	}

	// or even a hydronium!
	else if (num_H == 3) {
	  _mols.push_back (new Hydronium ());
	}

	else if (num_H == 0 || num_H > 3) {
	  printf ("XYZSystem::_ParseMolecules() - found water with %d H's\n", (int)atoms.size());
	  O->Print();
	  for (Atom_it jt = atoms.begin(); jt != atoms.end(); jt++) {
		(*jt)->Print();
	  }
	  exit(1);
	}

	MolPtr newmol = _mols[molIndex];
	newmol->MolID (molIndex);

	// let's set all the atom properties that we can, and add them into the molecule
	// and we also add in the hydrogens that are covalently bound - note, this can be more than 2 in the case of H3O+
	atoms.push_back(O);		// Don't forget to tack in the oxygen!
	for (Atom_it jt = atoms.begin(); jt != atoms.end(); jt++) {
	  newmol->AddAtom (*jt);
	}

	_UpdateUnparsedList(atoms);

	/*
	// now, from before, if one of the hydrogens is shared between two molecules making a contact-ion pair, then we will merge the two molecule to make one, and also update our ever-growing list of molecules to reflect it.
	if (s_mol != (Molecule *)NULL) {
	s_mol->Merge (newmol);		// this will swallow the new molecule into the contact ion pair
	delete newmol;				// get rid of the newmol
	vector<Molecule *>::iterator imol = _mols.end() - 1;	// fixes the _mols to get rid of the newmol we just made so it's not double-counted
	_mols.erase(imol);
	}
	 */
  }

  return;
}	// Parse Waters

void XYZSystem::_ParseNitrates () {

  /************************
   * Processing Nitric Acid
   * **********************/

  // The easiest thing to spot for nitric acid is, of course, the nitrogen! Only one of those, so let's grab it
  for (Atom_it N = _atoms.begin(); N != _atoms.end(); N++) {
	if ((*N)->Name().find("N") == string::npos) continue;

	// analogous to water, let's grab all the (covalently) bound O's of the molecule
	Atom_ptr_vec NAatoms (graph.BondedAtoms (*N, covalent, "O"));

	// at least 3 oxygens for a nitrate/nitric acid
	if (NAatoms.size() != 3) continue;
	bool fullNA = false;	// let's us know if the molecule is an hno3 or an no3
	AtomPtr no3H = (AtomPtr)NULL;				// this is the H bonded to an NO3 (we use it down below)
	AtomPtr no3O = (AtomPtr)NULL;				// this is the oxgen of the no3 that is closest to the H
	double no3OHdistance = 1000.0;				// the distance of the H to the nearest O

	// let's do a quick check to look at the oxygens of the molecule and see if it's a proper HNO3 or an NO3.
	for (Atom_it O = NAatoms.begin(); O != NAatoms.end(); O++) {

	  // these are all the Hs bound to the O
	  Atom_ptr_vec Hs (graph.BondedAtoms (*O, covalent, "H"));
	  if (!Hs.size()) continue;
	  // if an H is attached to one of the oxygens then we have a bonafide nitric acid
	  fullNA = true;

	  // here we'll run through all the Hs that are covalently bound and find which is the closest to an no3 oxygen.
	  for (Atom_it H = Hs.begin(); H != Hs.end(); H++) {
		double distance = graph.Distance (*O, *H);
		no3OHdistance = (distance < no3OHdistance) ? distance : no3OHdistance;
		if (distance == no3OHdistance) {
		  no3H = *H;
		  no3O = *O;
		}
	  }
	}

	int molIndex = (int)_mols.size();

	// and let's not forget to add the nitrogen
	NAatoms.push_back (*N);

	// if the molecule doesn't have a hydrogen - it's not an NA - it's a nitrate
	if (!fullNA) {
	  _mols.push_back (new Nitrate ());
	}

	// if it's a full nitric acid with proton:
	if (fullNA) {
	  NAatoms.push_back (no3H);
	  _mols.push_back (new NitricAcid());
	}

	_mols[molIndex]->MolID (molIndex);

	// let's set all the atom properties that we can, and add them into the molecule
	for (Atom_it O = NAatoms.begin(); O != NAatoms.end(); O++) {
	  _mols[molIndex]->AddAtom (*O);
	}

	_UpdateUnparsedList(NAatoms);

  }

  return;
}	// Parse Nitrates

void XYZSystem::_ParseSulfides () {
  for (Atom_it S = _atoms.begin(); S != _atoms.end(); S++) {
	if ((*S)->Name().find("S") == string::npos) continue;

	// for every S in the system, see if 2 oxygens are connected
	Atom_ptr_vec Os (graph.BondedAtoms (*S, covalent, "O"));
	if (Os.size() != 2) continue;

	int molIndex = (int)_mols.size();

	Os.push_back(*S);

	_mols.push_back (new SulfurDioxide());
	_mols[molIndex]->MolID (molIndex);

	for (Atom_it it = Os.begin(); it != Os.end(); it++) {
	  _mols[molIndex]->AddAtom (*it);
	}

	_UpdateUnparsedList(Os);
  }

}	// Parse SO2


void XYZSystem::_UpdateUnparsedList (const Atom_ptr_vec& parsed) {
  // fix up the tracking list
  for (Atom_ptr_vec::iterator it = _unparsed.begin(); it != _unparsed.end(); it++) {
	for (Atom_it jt = parsed.begin(); jt != parsed.end(); jt++) {
	  if (*it == *jt) {
		*it = (AtomPtr)NULL;
		break;
	  }
	}
  }
  return;
}
