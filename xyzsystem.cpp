#include "xyzsystem.h"

XYZSystem::XYZSystem (const std::string& filepath, const VecR& size, const std::string& wannierpath) :
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
	try {
		graph.UpdateGraph (_atoms);
	}
	catch (bondgraph::BondGraph::graphex& ex) {
		std::cout << "Caught an exception while updating the bond graph" << std::endl;
	}

	// Now let's do some house-cleaning to set us up for working with new molecules
	for (Mol_it it = _mols.begin(); it != _mols.end(); it++) {
		delete *it;				// get rid of all the molecules in memory
	}
	_mols.clear();				// then clear out the molecule list

	// we also have to go through and clear out some info on all the atoms
	// like the parentmolecules and names
	for (Atom_it it = _atoms.begin(); it != _atoms.end(); it++) {
		(*it)->ParentMolecule ( (MolPtr)NULL );
		(*it)->Residue ("");
		(*it)->MolID (-1);
	}

	// track which atoms have already been parsed in the system
	_unparsed.clear();
	std::copy (_atoms.begin(), _atoms.end(), std::back_inserter(_unparsed));

	this->_ParseSimpleMolecule<Hydroxide> (Atom::O, Atom::H, 1);
	this->_ParseSimpleMolecule<Water> (Atom::O, Atom::H, 2);
	this->_ParseSimpleMolecule<Hydronium> (Atom::O, Atom::H, 3);
	// parse NO3- ions
	this->_ParseSimpleMolecule<Nitrate> (Atom::N, Atom::O, 3);
	// turn any NO3- ions that have a covalently bound H into HNO3 - nitric acid
	this->_ParseNitricAcids ();
	this->_ParseSimpleMolecule<SulfurDioxide> (Atom::S, Atom::O, 2);
	this->_ParseProtons ();

	try {
		this->_CheckForUnparsedAtoms ();
	}
	catch (xyzsysex& ex) {
		std::cout << "Exception thrown while checking for unparsed atoms in the system - some atoms were not parsed into molecules" << std::endl;
		throw;
	}

	return;
}



void XYZSystem::_ParseNitricAcids () {

	for (Mol_it mol = _mols.begin(); mol != _mols.end(); mol++) {
		if ((*mol)->MolType() != Molecule::NO3) continue;

		AtomPtr N = (*mol)->GetAtom(Atom::N);
		// find any Hs bound to oxygens of the NO3
		Atom_ptr_vec Os = graph.BondedAtoms(N, bondgraph::covalent, Atom::O);
		Atom_ptr_vec Hs; 
		for (Atom_it O = Os.begin(); O != Os.end(); O++) {
			Atom_ptr_vec H_check = graph.BondedAtoms(*O, bondgraph::covalent, Atom::H);
			if (H_check.size() == 1)
				Hs.push_back(H_check[0]);
			if (H_check.size() > 1) {
				printf ("A nitrate Oxygen has %zu covalently bound hydrogens! what's going on here?\n", H_check.size());
				(*O)->Print();
				exit(1);
			}
		}

		if (Hs.size() == 1)
			(*mol)->AddHydrogen (Hs[0]);
		else if (Hs.size() > 1) {
			printf ("This nitric acid has %zu covalently bound hydrogens!\n", Hs.size());
			(*mol)->Print();
			for (Atom_it ai = Hs.begin(); ai != Hs.end(); ai++) {
				(*ai)->Print();
			}
		}
	}
	return;
}	// Parse Nitric acids


void XYZSystem::_UpdateUnparsedList (Atom_ptr_vec& parsed) {
	// fix up the list for tracking atoms that have already been added into molecules.
	// _unparsed should contain only atoms that are not in the parsed list
	// these atom vectors get sorted according the the atom's ID

	_unparsed.erase(
			std::remove_if (_unparsed.begin(), _unparsed.end(), std::bind2nd(AtomPtr_In_List<Atom_ptr_vec>(), parsed)), _unparsed.end());

	return;
}


void XYZSystem::_CheckForUnparsedAtoms () const {

	if (!_unparsed.empty()) {

		std::cout << "The following atoms were found unparsed into molecules after all molecules had been formed" << std::endl;

		// print out every atom that hasn't been parsed
		for (Atom_it it = _unparsed.begin(); it != _unparsed.end(); it++) {
			std::cout << std::endl;
			(*it)->Print();	// show all remaining atoms

			// and all the atoms to which it is bound
			Atom_ptr_vec bound (graph.BondedAtoms(*it));
			for (Atom_it jt = bound.begin(); jt != bound.end(); jt++) { 
				// and the distance between them
				std::cout << "^--~ (" << graph.Distance(*jt, *it) << ")  ";
				(*jt)->Print();
			}
		}
		throw (unaccountedex());
	}

	return;
}	// check for unparsed atoms




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

	for (Atom_it it = _atoms.begin(); it != _atoms.end(); it++) {
		// find every oxygen and sulfur atom (those are the ones that contain wanniers
		if ((*it)->Element() != Atom::O && (*it)->Element() != Atom::S) continue;

		MolPtr mol = (*it)->ParentMolecule();

		// and then go through all the wanniers in order to find the ones attached to the oxygen
		for (VecR_it vi = _wanniers.begin(); vi != _wanniers.end(); vi++) {

			// we find the distance to the oxygen from the wannier center
			double distance = MDSystem::Distance ((*it)->Position(), *vi).Magnitude();

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

	try {
		this->_ParseMolecules();
	} catch (xyzsysex& ex) {
		std::cout << "Exception caught while parsing the molecules of the XYZ system" << std::endl;
		throw;
	}

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

	try {
		if (++_reparse_step == _reparse_limit) {
			this->_ParseMolecules();
			_reparse_step = 0;
		}
	} catch (xyzsysex& ex) {
		std::cout << "Exception caught while parsing the molecules of the XYZ system" << std::endl;
		throw;
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

	for (Atom_it it = _atoms.begin(); it != _atoms.end(); it++) {
		VecR ri = (*it)->Position();
		dipole += ri * (*it)->Charge();
	}

	for (VecR_it it = _wanniers.begin(); it != _wanniers.end(); it++) {
		dipole -= (*it) * 2.0;
	}

	return (dipole);
}



void XYZSystem::_ParseProtons () {

	Atom_ptr_vec parsed;

	for (Atom_it H = _atoms.begin(); H != _atoms.end(); H++) {
		// form all the remaining un-bound hydrogens into their own "proton" molecules
		if ((*H)->Element() != Atom::H) continue;

		Atom_ptr_vec cov = graph.BondedAtoms(*H, bondgraph::covalent);
		if (cov.empty()) {

			int molIndex = (int)_mols.size();	// set the molecule index
			_mols.push_back (new Proton ());

			MolPtr newmol = _mols[molIndex];
			newmol->MolID (molIndex);
			newmol->AddAtom (*H);

			parsed.push_back (*H);
		}
		// otherwise, if there are covalent bonds, AND the hydrogen is not parsed into a molecule... something is wrong
		else if (!(*H)->ParentMolecule()) {
			printf ("Found a hydrogen that is not parsed into a molecule, but has covalent bonds!\n");
			(*H)->Print();
			printf ("It is covalently bound to:\n");
			for (Atom_it it = cov.begin(); it != cov.end(); it++) {
				(*it)->Print();
			}
			printf ("And it is h-bonded to:\n");
			Atom_ptr_vec hbonds = graph.BondedAtoms(*H, bondgraph::hbond);
			for (Atom_it it = hbonds.begin(); it != hbonds.end(); it++) {
				(*it)->Print();
			}

			exit(1);
		}
	}

	_UpdateUnparsedList (parsed);
	return;
}



