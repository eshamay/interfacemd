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
	_bondgraph.UpdateGraph (_atoms.Atoms());

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
	
/*******************
 * Processing Waters
 *******************/

	// let's go through the system and find all the waters and form molecules out of them
	RUN (_atoms) {

		// every water has an oxygen atom, right? So let's find those
		Atom * O = _atoms[i];
		if (O->Name().find("O") == string::npos) continue;

		// The bondgraph provides us with all the covalently bound H's
		vector<Atom *> atoms = _bondgraph.CovalentBonds (O);
		vector<Atom *> Hs = _bondgraph.CovalentBonds (O, "H");

		// if we pick up an OH group that is part of a larger molecule (i.e. nitric acid) then it will be processed as a hydroxide...
		// So if the only type of atom attached is a hydrogen, we have some form of water
		if (atoms.size() != Hs.size()) continue;

		int molIndex = _mols.size();	// set the molecule index

		// if any of the Hs are being shared between two molecules (i.e. a contact-ion pair) then we have to resolve which molecule gets the atom.
		// A simple solution is to decide that the molecule with the oxygen closer to the hydrogen is the one that wins out
		atoms.clear();	// we'll add in all the non-contentious hydrogens here

		RUN (Hs) {
			Atom * H = Hs[i];
			// if the hydrogen hasn't yet been assigned, then there is no contention - add it into the current molecule
			if (H->ParentMolecule() == (Molecule *)NULL) {
				atoms.push_back(H);
			}
			// otherwise, we have to find the distance from the current oxygen to the hydrogen, and the other oxygen to the hydrogen
			else {
				Molecule * wat2 = H->ParentMolecule();
				Atom * O2 = wat2->GetAtom("O");
				double distance1 = _bondgraph.Distance (H, O);
				double distance2 = _bondgraph.Distance (H, O2);
				
				// if the H is closer to the current oxygen then we remove it from the other molecule, and add it into this one
				if (distance1 < distance2) {
					wat2->RemoveAtom (H);
					atoms.push_back(H);
				}
				// otherwise... we just jump over this hydrogen and leave it in the other molecule
			}
		}

	  	// we may be dealing with a hydroxide ion
	  	if (atoms.size() == 1) {
	  		_mols.push_back (new Hydroxide ());
	  	}

	  	// otherwise we have a full molecule
	  	if (atoms.size() == 2) {
			_mols.push_back (new Water ());
	  	}

		// or even a hydronium!
		if (atoms.size() == 3) {
			_mols.push_back (new Hydronium ());
		}

		if (atoms.size() == 0) { 
			//printf ("found water with %d H's\n", Hcount); 
			continue; 
		}

		_mols[molIndex]->MolID (molIndex);
			

		// let's set all the atom properties that we can, and add them into the molecule
		// and we also add in the hydrogens that are covalently bound - note, this can be more than 2 in the case of H3O+
		atoms.push_back(O);		// lastly, tack in the oxygen!
		RUN (atoms) {
			_mols[molIndex]->AddAtom (atoms[i]);
		}
	}

/************************
 * Processing Nitric Acid
 * **********************/

	// The easiest thing to spot for nitric acid is, of course, the nitrogen! Only one of those, so let's grab it
	for (int N = 0; N < _atoms.size(); N++) {
		if (_atoms[N]->Name().find("N") == string::npos) continue;

		// analogous to water, let's grab all the (covalently) bound O's of the molecule
		std::vector<Atom *> NAatoms = _bondgraph.CovalentBonds (_atoms[N], "O");
		
		int Ocount = NAatoms.size();	// number of Os on the molecule
		// at least 3 oxygens for a nitrate/nitric acid
		if (Ocount != 3) continue;

		bool fullNA = false;	// let's us know if the molecule is an hno3 or an no3
		Atom * no3H = (Atom *)NULL;				// this is the H bonded to an NO3 (we use it down below)
		Atom * no3O = (Atom *)NULL;				// this is the oxgen of the no3 that is closest to the H 
		double no3OHdistance = 1000.0;				// the distance of the H to the nearest O

		// let's do a quick check to look at the oxygens of the molecule and see if it's a proper HNO3 or an NO3.
		RUN (NAatoms) {
			
			Atom * O = NAatoms[i];

			// these are all the Hs bound to the O
			std::vector<Atom *> Hs (_bondgraph.CovalentBonds (O, "H"));

			if (!Hs.size()) continue;

			// here we'll run through all the Hs that are covalently bound and find which is the closest to an no3 oxygen.
			RUN (Hs) {
				double distance = _bondgraph.Distance (O, Hs[i]);
				no3OHdistance = (distance < no3OHdistance) ? distance : no3OHdistance;
				fullNA = true;
				if (distance == no3OHdistance) {
					no3H = Hs[i];
					no3O = O;
				}
			}
		}

	
		// now... there are situations where the no3 and its closest water are actually sharing the hydrogen. We have to check here to see if it is closer to the water or closer to the no3 group to assign it either way.
		// first we check if there was even an H in the molecule
		if (fullNA) {
			// if there is an H, let's check to see if it's already assigned to a water from before:
	
			// first, is the H even attached to another molecule?
			if (no3H->ParentMolecule() != (Molecule *)NULL) {
				// if it is shared, then we have to ***ASSUME*** that it's being shared with a water molecule (for now).
				// find the distance to the bonded oxygen of the water
				Water * wat = static_cast<Water *>(no3H->ParentMolecule());
				Atom * watO = wat->GetAtom("O");
				double watOHdistance = _bondgraph.Distance (no3H, watO);
		
				// If the H is closer to the no3, then we call this a full hno3, and move the atom from the water to this new molecule.
				if (no3OHdistance < watOHdistance) {
					// first remove it from the water
					wat->RemoveAtom (no3H);
					// then add it in to be processed into the new nitric acid
					fullNA = true;
				}
				// otherwise, everything stays the same, and we treat the new molecule as a nitrate
				else {
					fullNA = false;
				}
			}
			// if the H is not contended between two molecules, then the no3 gets it!
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
		_wanniers.LoadFirst(); 
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
