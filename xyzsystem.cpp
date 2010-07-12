#include "xyzsystem.h"

XYZSystem::XYZSystem (string filepath, VecR size, string wannierpath) :
	MDSystem (),
	_coords(filepath),
	_wanniers(wannierpath)
{

	MDSystem::Dimensions (size);

	this->LoadFirst();
}

XYZSystem::~XYZSystem () {
	RUN (_mols) {
		delete _mols[i];
	}
}

void XYZSystem::_ParseMolecules () {

/***********************************************************************************
 * This is the top-level parsing routine to give the overall idea of what's going on
 * *********************************************************************************/

	// first things first - we need the interatomic distances!
	_graph.UpdateGraph (_atoms);

	// Now let's do some house-cleaning to set us up for working with new molecules
	RUN (_mols)
		delete _mols[i];		// get rid of all the molecules in memory
	_mols.clear();				// then clear out the molecule list

	// we also have to go through and clear out some info on all the atoms
	// like the parentmolecules and names
	RUN (_atoms) {
		_atoms[i]->ParentMolecule ( (Molecule *)NULL );
		_atoms[i]->Residue ("");
		_atoms[i]->MolID (-1);
	}

	// track which atoms have already been parsed in the system
	_unparsed.clear();
	_unparsed = _atoms;

	this->_ParseWaters ();
	this->_ParseNitrates ();

	bool leave = false;
	RUN (_unparsed) {
		Atom * a1 = _unparsed[i];
		if (a1 != (Atom *)NULL) {
			a1->Print();		// show all remaining atoms
			std::vector<Atom *> bound (_graph.BondedAtoms(a1));
			RUN2(bound) {		// and all the atoms to which it is bound (and the distance)
				cout << "^--~ (" << _graph.Distance(bound[j], a1) << ")  ";
				bound[j]->Print();
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
	RUN (_mols) {
		_mols[i]->ClearWanniers();
	}

	// try a new approach here
	vector<VecR> wans (_wanniers.Coords());
	vector<VecR>::iterator vi;

	RUN (_atoms) {
		// find every oxygen atom
		if (_atoms[i]->Name().find("O") == string::npos) continue;

		Atom * O = _atoms[i];
		Molecule * mol = O->ParentMolecule();

		// only load up water wanniers for now...
		//if (mol->Name().find("h2o") == std::string::npos) continue;

		// and then go through all the wanniers in order to find the ones attached to the oxygen
		for (vi = wans.begin(); vi != wans.end(); vi++) {

			// we find the distance to the oxygen from the wannier center
			double distance = MDSystem::Distance (O->Position(), *vi).Magnitude();

			// if it's close enough, then that oxygen gets it
			if (distance < 1.0) {
				mol->AddWannier(*vi);
				// and we take it out of the list so that we only scan the remaining wanniers for the remaining atoms and save a little time
			}

		}
	}

/*
	//bool home = false;

	// Right now we calculate the wannier centers for each molecule based on the oxygens. So run through each oxygen and find the 4 nearest wanniers
	RUN (_atoms) {
		if (_atoms[i]->Name().find("O") == string::npos) continue;

		vector< vector<double> > distances;

		// run through and find the distance to each wannier center
		for (int wan = 0; wan < _wanniers.size(); wan++) {
			vector<double> temp;

			double distance = _atoms[i]->Position().MinDistance (_wanniers[wan], Atom::Size());
			temp.push_back(distance);
			temp.push_back((double)wan);

			distances.push_back(temp);
		}

		// then sort the wannier centers by distance to the oxygen
		sort(distances.begin(), distances.end());

		//printf ("setting wan in oxy (%d)\t parent = %s  [%d]\n",
			//_atoms[i]->ID(), _atoms[i]->Residue().c_str(), _atoms[i]->MolID());
		// then add in the 4 closest wannier centers
		for (int wan = 0; wan < 4; wan++) {
			int id = int(distances[wan][1]);
			_atoms[i]->ParentMolecule()->AddWannier(_wanniers[id]);
		}
	}
*/
return;
}

void XYZSystem::LoadFirst () {
	//_atoms.LoadFirst();
	_atoms.clear();
	_atoms = _coords.Atoms();
	this->_ParseMolecules();

	if (_wanniers.Loaded()) {
		//_wanniers.LoadFirst();
		this->_ParseWanniers();
	}

}

// watch out - no functionality for wannier centers here (yet)
void XYZSystem::Seek (int step) {
	_coords.Seek (step);
	_atoms.clear();
	_atoms = _coords.Atoms();
	this->_ParseMolecules();
}

void XYZSystem::LoadNext () {
	_coords.LoadNext();
	_atoms.clear();
	_atoms = _coords.Atoms();
	this->_ParseMolecules();

	if (_wanniers.Loaded()) {
		_wanniers.LoadNext();
		this->_ParseWanniers();
	}
}

// This will calculate the total dipole moment of the system based on atom locations and wannier center positions
// The origin is shifted to the center of the system in order to get closest images (wrapped into the box) of all the atoms/wanniers
VecR XYZSystem::SystemDipole () {

	VecR dipole;
	//VecR center = Atom::Size() * 0.5;	// center of the system

	RUN (_atoms) {
		//VecR ri = center.MinVector(_atoms[i]->Position(), Atom::Size());
		VecR ri = _atoms[i]->Position();
		dipole += ri * _atoms[i]->Charge();
	}

	RUN (_wanniers) {
		//VecR ri = center.MinVector(_wanniers[i], Atom::Size());
		VecR ri = _wanniers[i];
		dipole -= ri * 2.0;
	}

return (dipole);
}

void XYZSystem::_ParseWaters () {

/*******************
 * Processing Waters
 *******************/

	// let's go through the system and find all the waters and form molecules out of them
	RUN (_atoms) {

		// every water has an oxygen atom, right? So let's find those
		Atom * O = _atoms[i];
		if (O->Name().find("O") == string::npos) continue;

		// The bondgraph provides us with all the covalently bound H's
		vector<Atom *> atoms = _graph.BondedAtoms (O, covalent);
		vector<Atom *> Hs = _graph.BondedAtoms (O, covalent, "H");

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
			RUN (atoms) {
				atoms[i]->Print();
			}
			exit(1);
		}

		Molecule * newmol = _mols[molIndex];
		newmol->MolID (molIndex);

		// let's set all the atom properties that we can, and add them into the molecule
		// and we also add in the hydrogens that are covalently bound - note, this can be more than 2 in the case of H3O+
		atoms.push_back(O);		// Don't forget to tack in the oxygen!
		RUN (atoms) {
			newmol->AddAtom (atoms[i]);
		}

#ifdef DEBUG
		// fix up the tracking list
		RUN (_unparsed) {
			Atom * a1 = _unparsed[i];
			RUN2 (atoms) {
				if (a1 == atoms[j]) {
					_unparsed[i] = (Atom *)NULL;
					break;
				}
			}
		}
#endif

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
}

void XYZSystem::_ParseNitrates () {

/************************
 * Processing Nitric Acid
 * **********************/

	// The easiest thing to spot for nitric acid is, of course, the nitrogen! Only one of those, so let's grab it
	for (size_t N = 0; N < _atoms.size(); N++) {
		if (_atoms[N]->Name().find("N") == string::npos) continue;

		// analogous to water, let's grab all the (covalently) bound O's of the molecule
		std::vector<Atom *> NAatoms (_graph.BondedAtoms (_atoms[N], nobond, "O"));

		// at least 3 oxygens for a nitrate/nitric acid
		if (NAatoms.size() != 3) continue;
		bool fullNA = false;	// let's us know if the molecule is an hno3 or an no3
		Atom * no3H = (Atom *)NULL;				// this is the H bonded to an NO3 (we use it down below)
		Atom * no3O = (Atom *)NULL;				// this is the oxgen of the no3 that is closest to the H
		double no3OHdistance = 1000.0;				// the distance of the H to the nearest O

		// let's do a quick check to look at the oxygens of the molecule and see if it's a proper HNO3 or an NO3.
		RUN (NAatoms) {

			Atom * O = NAatoms[i];

			// these are all the Hs bound to the O
			std::vector<Atom *> Hs (_graph.BondedAtoms (O, covalent, "H"));
			if (!Hs.size()) continue;
			// if an H is attached to one of the oxygens then we have a bonafide nitric acid
			fullNA = true;

			// here we'll run through all the Hs that are covalently bound and find which is the closest to an no3 oxygen.
			RUN (Hs) {
				double distance = _graph.Distance (O, Hs[i]);
				no3OHdistance = (distance < no3OHdistance) ? distance : no3OHdistance;
				if (distance == no3OHdistance) {
					no3H = Hs[i];
					no3O = O;
				}
			}
		}

		int molIndex = (int)_mols.size();

		// and let's not forget to add the nitrogen
		NAatoms.push_back (_atoms[N]);

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
		RUN (NAatoms) {
			_mols[molIndex]->AddAtom (NAatoms[i]);
		}

		// fix up the tracking list
		RUN (_unparsed) {
			Atom * a1 = _unparsed[i];
			RUN2 (NAatoms) {
				if (a1 == NAatoms[j]) {
					_unparsed[i] = (Atom *)NULL;
					break;
				}
			}
		}

	}

return;
}
