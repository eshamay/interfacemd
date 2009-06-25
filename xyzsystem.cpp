#include "xyzsystem.h"

XYZSystem::XYZSystem (string filepath, VecR size, string wannierpath)
	: _atoms(XYZFile(filepath)), _dims(size), _wanniers(WannierFile(wannierpath))
{

	// set the system size
	Atom::Size (size);

	this->LoadFirst();
}

XYZSystem::~XYZSystem () {

	RUN (_mols) {
		delete _mols[i];
	}

	RUN (_atoms) {
		delete _atoms[i];
	}
}

void XYZSystem::_ParseMolecules () {

/***********************************************************************************
 * This is the top-level parsing routine to give the overall idea of what's going on
 * *********************************************************************************/

	// first things first - we need the interatomic distances!
	_matrix.UpdateMatrix (_atoms.Atoms());

	// Now let's do some house-cleaning to set us up for working with new molecules
	RUN (_mols) {
		delete _mols[i];		// get rid of all the molecules in memory
	}
	_mols.clear();				// then clear out the molecule list

	// we also have to go through and clear out some info on all the atoms
	// like the parentmolecules and names
	RUN (_atoms) {
		_atoms[i]->ParentMolecule ( (Molecule *)NULL );
		_atoms[i]->Residue ("");
		_atoms[i]->MolID (-1);
	}

	/* For debugging (and other useful things?) this will keep a list of all the atoms that have been processed into molecules. Any atoms left over at the end of the parsing routine are not included and ... can potentially cause problems */
	std::vector<Atom *> atom_tracking_list (_atoms.Atoms());

/*******************
 * Processing Waters
 *******************/

	// let's go through the system and find all the waters and form molecules out of them
	RUN (_atoms) {

		// every water has an oxygen atom, right? So let's find those
		Atom * O = _atoms[i];
		if (O->Name().find("O") == string::npos) continue;

		// The bondgraph provides us with all the covalently bound H's
		vector<Atom *> atoms = _matrix.BondedAtoms (O, covalent);
		vector<Atom *> Hs = _matrix.BondedAtoms (O, covalent, "H");

		// if we pick up an OH group that is part of a larger molecule (i.e. nitric acid) then it will be processed as a hydroxide...
		// So if the only type of atom attached is a hydrogen, we have some form of water (OH, H2O, H3O)
		if (Hs.size() != atoms.size()) continue;

		int molIndex = _mols.size();	// set the molecule index

		// if any of the Hs are being shared between two molecules (i.e. a contact-ion pair) then we have to resolve which molecule gets the atom.
		// A simple solution is to decide that the molecule with the oxygen closer to the hydrogen is the one that wins out
		atoms.clear();	// we'll add in all the non-contentious hydrogens here

		Molecule * s_mol = (Molecule *)NULL;
		RUN (Hs) {
			Atom * H = Hs[i];
			s_mol = H->ParentMolecule();
			// if the hydrogen hasn't yet been assigned, then there is no contention - add it into the current molecule
			if (s_mol == (Molecule *)NULL) {		// atom is still available - not connected
				atoms.push_back(H);
			}
			// otherwise we have to find the molecule that it's closer to, and remove it from the other
			else {
				Atom * s_O = s_mol->GetAtom("O");		// **assuming** that it's bound to an O
				double t_distance = _matrix.Distance(H, O);		// distance to target O ("this" O)
				double s_distance = _matrix.Distance(H, s_O);	// distance to source O ("other" O)
				// now, if the distance between the H and the other molecule's atom is smaller, then we change a few things. Namely, we have to set the bond between the H and this current O to hydrogen, instead of covalent...
				if (s_distance < t_distance) {
					//_matrix.SetBond(H, O, hbond);
				#ifdef DEBUG
					cout << "XYZSystem::_ParseMolecules() - Found a shared H!" << endl;
					H->Print();
					cout << "it's closer to:" << endl;
					s_mol->Print();
					cout << "and shared with:" << endl;
					O->Print();
				#endif
				}
				// otherwise, if the distance to the target O (the current one) is smaller, then we have to set the other bond to an H-bond, remove the H from that molecule, and add it into the current one
				else {
					//_matrix.SetBond(H, s_O, hbond);
					s_mol->RemoveAtom (H);
					atoms.push_back(H);

				#ifdef DEBUG
					cout << "XYZSystem::_ParseMolecules() - Found a shared H!" << endl;
					H->Print();
					cout << "it's closer to:" << endl;
					O->Print();
					cout << "and shared with:" << endl;
					s_mol->Print();
				#endif
				}
			}
		}

	  	// we may be dealing with a hydroxide ion
	  	if (atoms.size() == 1) {
	  		_mols.push_back (new Hydroxide ());
	  	}

	  	// otherwise we have a full molecule
	  	else if (atoms.size() == 2) {
			_mols.push_back (new Water ());
	  	}

		// or even a hydronium!
		else if (atoms.size() == 3) {
			_mols.push_back (new Hydronium ());
		}

		else if (atoms.size() == 0) {
			printf ("XYZSystem::_ParseMolecules() - found water with %d H's\n", atoms.size());
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

		// fix up the tracking list
		RUN (atom_tracking_list) {
			Atom * a1 = atom_tracking_list[i];
			if (a1 == O) {
				atom_tracking_list[i] = (Atom *)NULL;
				continue;
			}
			RUN2 (atoms) {
				if (a1 == atoms[j]) {
					atom_tracking_list[i] = (Atom *)NULL;
					break;
				}
			}
		}

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

/************************
 * Processing Nitric Acid
 * **********************/

	// The easiest thing to spot for nitric acid is, of course, the nitrogen! Only one of those, so let's grab it
	for (int N = 0; N < _atoms.size(); N++) {
		if (_atoms[N]->Name().find("N") == string::npos) continue;

		// analogous to water, let's grab all the (covalently) bound O's of the molecule
		std::vector<Atom *> NAatoms (_matrix.BondedAtoms (_atoms[N], nobond, "O"));

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
			std::vector<Atom *> Hs (_matrix.BondedAtoms (O, covalent, "H"));
			if (!Hs.size()) continue;
			// if an H is attached to one of the oxygens then we have a bonafide nitric acid
			fullNA = true;

			// here we'll run through all the Hs that are covalently bound and find which is the closest to an no3 oxygen.
			RUN (Hs) {
				double distance = _matrix.Distance (O, Hs[i]);
				no3OHdistance = (distance < no3OHdistance) ? distance : no3OHdistance;
				if (distance == no3OHdistance) {
					no3H = Hs[i];
					no3O = O;
				}
			}
		}


		// now... there are situations where the no3 and its closest water are actually sharing the hydrogen. We have to check here to see if it is closer to the water or closer to the no3 group to assign it either way.
		if (fullNA) {
			// if there is an H, let's check to see if it's already assigned to another molecule from before:
			Molecule * s_mol = no3H->ParentMolecule();

			// first, is the H even attached to another molecule?
			if (s_mol != (Molecule *)NULL) {
				#ifdef DEBUG
					cout << "Parsing an NO3 and Found an H shared with:" << endl;
					s_mol->Print();
				#endif
				// if it is shared, then we have to find which it is closer to.
				// here we ***assume*** that the H is bound to an oxygen
				Atom * s_O = s_mol->GetAtom("O");
				double s_distance = _matrix.Distance (no3H, s_O);

				// If the H is closer to the no3, then:
				if (no3OHdistance < s_distance) {
				// 		1) fix the source bond to be an H-bond instead of covalent
					//_matrix.SetBond(s_O, no3H, hbond);
				// 		2) remove the H from the other molecule
					s_mol->RemoveAtom (no3H);
				// 		3) add the H into this one
					fullNA = true;
				}
				// otherwise :
				else {
				// 		1) fix the bondtype in the adjacency matrix to make it an H-bond
					//_matrix.SetBond(no3H, no3O, hbond);
				// 		2) turn this current molecule into a nitrate
					fullNA = false;
				}
			}
			// if the H is not contended between two molecules, then the "current" no3 gets it!
			else {
				fullNA = true;
			}
		}

		unsigned int molIndex = _mols.size();

		// if the molecule doesn't have a hydrogen - it's not an NA - it's a nitrate
		if (!fullNA) {
			_mols.push_back (new Nitrate ());
		}

		// if it's a full nitric acid with proton:
		if (fullNA) {
			_mols.push_back (new NitricAcid());
			NAatoms.push_back (no3H);
		}

		// and let's not forget to add the nitrogen
		NAatoms.push_back (_atoms[N]);

		_mols[molIndex]->MolID (molIndex);

		// let's set all the atom properties that we can, and add them into the molecule
		RUN (NAatoms) {
			_mols[molIndex]->AddAtom (NAatoms[i]);
		}

		// fix up the tracking list
		RUN (atom_tracking_list) {
			Atom * a1 = atom_tracking_list[i];
			RUN2 (NAatoms) {
				if (a1 == NAatoms[j]) {
					atom_tracking_list[i] = (Atom *)NULL;
					break;
				}
			}
		}

	}

	bool leave = false;
	RUN (atom_tracking_list) {
		Atom * a1 = atom_tracking_list[i];
		if (a1 != (Atom *)NULL) {
			a1->Print();		// show all remaining atoms
		#ifdef DEBUG
			std::vector<Atom *> bound (_matrix.BondedAtoms(a1));
			RUN2(bound) {		// and all the atoms to which it is bound (and the distance)
				cout << "^--~ (" << _matrix.Distance(bound[j], a1) << ")  ";
				bound[j]->Print();
			}
		#endif
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
			double distance = O->Position().MinDistance (*vi, Atom::Size());

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
	_atoms.LoadFirst();

	this->_ParseMolecules();

	if (_wanniers.Loaded()) {
		//_wanniers.LoadFirst();
		this->_ParseWanniers();
	}

}

// watch out - no functionality for wannier centers here (yet)
void XYZSystem::Seek (int step) {
	_atoms.Seek (step);
	this->_ParseMolecules();
}

void XYZSystem::LoadNext () {
	_atoms.LoadNext();
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
	VecR center = Atom::Size() * 0.5;	// center of the system

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
