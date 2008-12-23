#include "ambersystem.h"

AmberSystem::AmberSystem (string prmtop, string mdcrd, string mdvel = "") 
	// some initialization needs to happen here
	: 	_topfile(TOPFile(prmtop)), 
		_coords(CRDFile(mdcrd, _topfile.NumAtoms())), 
		_forces(ForceFile (mdvel, _topfile.NumAtoms()))
{

	// because some really useful functionality comes out of the Atom class if the Atom::Size() is set, we'll do that here
	Atom::Size (_coords.Dims());

	// Then we have to form proper atoms. First setup the vector to hold them all.
	_atoms.resize(_topfile.NumAtoms()); // = vector<Atom> (_topfile.NumAtoms(), Atom());
	// and then actually create these bad mamma jammas
	RUN (_atoms) {
		_atoms[i] = new Atom ();
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
}

// While the crdfile holds spatial coordinate information, and the topology file holds atomic information, the data has to be processed into proper atoms in order to play around with them more effectively.
// This function should only be used once when first loading the system up. This info doesn't change during the course of the MD run
void AmberSystem::_ParseAtomInformation () {
	RUN (_atoms) {
		_atoms[i]->Name (_topfile.AtomNames()[i]);
		_atoms[i]->SetMass();
		_atoms[i]->SetCharge();
		_atoms[i]->ID (i);		// set the atom's index number - because we may need to access ordered/list info elsewhere
	}

return;
}

// This function copies the coordinates and force vectors into each of the atoms from the input files
void AmberSystem::_ParseAtomVectors () {

	RUN (_atoms) {
		_atoms[i]->Position ( _coords[i] );
		if (_forces.Loaded()) _atoms[i]->Force ( _forces[i] );
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
		string name = _topfile.MolNames()[mol];
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
		// otherwise add on a generic molecule
		else _mols.push_back (new Molecule());

		_mols[mol]->Name (name);	// set the molecule's residue name

		// Then add all the atoms between the indexes of the molpointers in the topology file
		// MPI had a weird problem **HERE**... fixed June 2007 ~ESS
		for (int atomCount = 0; atomCount < _topfile.MolSizes()[mol]; atomCount++) {
			// sets the index of the current atom that's being added to the molecule
			int curAtom = _topfile.MolPointers()[mol] + atomCount - 1;	

			// Now we're going to bless this new atom with loads of information about itself and its molecule
			_mols[mol]->AddAtom( _atoms[curAtom] );			// First add the atom into the molecule
			_mols[mol]->MolID (mol);
			_atoms[curAtom]->Residue (_mols[mol]->Name());	// While we're here, let's set the residue names
			_atoms[curAtom]->MolID (mol);					// set the molecule's ID #
			_atoms[curAtom]->ParentMolecule (_mols[mol]);	// Hell, why not also let the atom know what molecule it's in!
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
	RUN (_atoms)
		_atoms[i]->ParentMolecule()->Unset ();
	
return;
}

void AmberSystem::PrintCRDFile (string filepath) {

	RUN (_atoms) {
		printf ("  % 4.7f  % 4.7f  % 4.7f", _atoms[i]->Position()[x], _atoms[i]->Position()[y], _atoms[i]->Position()[z]);
	}
	printf ("  % 4.7f  % 4.7f  % 4.7f", this->Dims()[x], this->Dims()[y], this->Dims()[z]);

return;
}
