#include "connectmatrix.h"

const string CoordName (coordination coord) { 
	return COORD_NAMES[coord]; 
}

ConnectivityMatrix::ConnectivityMatrix (std::vector<Atom *>& atoms) : _atoms(atoms) { 

	// and resize the matrix to include all the atoms
	_matrix.resize(_atoms.size());

	RUN (_matrix) {
		_matrix[i].resize(_atoms.size());
		_atoms[i]->ClearHBonds();
	}

	this->UpdateMatrix();
}

// this fills in the bottom triangle of the matrix. A value of 1.0 means we have a covalent bond
void ConnectivityMatrix::_FormCovalentBond (const Atom * atom1, const Atom * atom2) {

	int index1 = atom1->ID();
	int index2 = atom2->ID();

	if (index1 > index2) _matrix[index1][index2] = 1.0;
	else _matrix[index2][index1] = 1.0;

return;
}

// let's us update the H-bond information (diagonal elements) and also let the atoms know they are H-bound to something
void ConnectivityMatrix::_FormHBond (Atom * atom1, Atom * atom2) {	

	_matrix[atom1->ID()][atom1->ID()]++;
	_matrix[atom2->ID()][atom2->ID()]++;

	atom1->FormHBond (atom2);
	atom2->FormHBond (atom1);

return;
}

// here we form a list of all the covalent bonds an atom is involved in
std::vector<Atom *> ConnectivityMatrix::CovalentBonds (const Atom * atom) const {

	std::vector<Atom *> output;
	int index = atom->ID();

	// first check the row
	for (int i = 0; i < index; i++) {
		if (_matrix[index][i] == 1.0) output.push_back (_atoms[i]);
	}

	// then check the column
	for (unsigned int i = index+1; i < _atoms.size(); i++) {
		if (_matrix[i][index] == 1.0) output.push_back (_atoms[i]);
	}

	return output;
}

// returns the distance between two atoms. Note that the first index used has to be the smaller of the two if we're going to access the top half of the matrix
double ConnectivityMatrix::Distance (const Atom * atom1, const Atom * atom2) const {

	int index1 = atom1->ID();
	int index2 = atom2->ID();

	double distance = (index1 < index2) ? _matrix[index1][index2] : _matrix[index2][index1];

return(distance);
}

// goes through and find all the bonds in the system based on the criteria above
void ConnectivityMatrix::_FindBonds (Atom * atom1, Atom * atom2) {


	double distance = this->Distance(atom1, atom2);
	//if (distance < H_BOND_LENGTH && distance > OH_BOND_LENGTH) printf ("%d\t%d\n", atom1->ID(), atom2->ID());

	// bonding between O and H (if either combo of H1 and O2 or H2 and O1 - we have an O...H pair)
	if ( (atom1->Name().find("H") != string::npos && atom2->Name().find("O") != string::npos)
		|| (atom2->Name().find("H") != string::npos && atom1->Name().find("O") != string::npos) ) {

		// let's form placeholders 
		Atom * o, * h;
		//Atom * o2;
		if (atom1->Name().find("H") != string::npos) { h = atom1; o = atom2; }
		else { h = atom2; o = atom1; }
//		o2 = (*h->ParentMolecule())["O"];
			
		// check here the covalent OH-oscillator criteria (distance)
		if (distance < OH_BOND_LENGTH) 
			this->_FormCovalentBond(o, h);
		
		// check the H-bond criteria (larger than a covalent OH-bond but smaller than the H-bond cutoff)
		if (distance < H_BOND_LENGTH && distance > OH_BOND_LENGTH) {
			this->_FormHBond (o, h);
		
/* 2ndly, check the angle cut-off criteria. This is explained well in the paper by DSW & Dennis Hore (J Phys Chem B, vol 110, No 41, 2006). To summarize: an H-bond is composed of an oxygen and hydrogen that are not covalently bound. Of that pair, find the covalent OH vector of the hydrogen, and the vector of the h-bonding pair. The angle between those two vectors should be < 30 degrees (according to the paper). */
/*
			o2 = (*h->ParentMolecule())["O"];
			// this is the covalent oh bond
			VecR oh = o2->Position().MinVector(h->Position(), Atom::Size());
			// this is the hydrogen bond
			VecR o_h = h->Position().MinVector(o->Position(), Atom::Size());

			if ((oh < o_h) < H_BOND_ANGLE) {
				this->_FormHBond (o, h);
			}
*/
		}
	}
		
	// bonding between N and O
	if ( (atom1->Name().find("N") != string::npos && atom2->Name().find("O") != string::npos)
		|| (atom2->Name().find("N") != string::npos && atom1->Name().find("O") != string::npos) ) {

		if (distance < NO_BOND_LENGTH) 
			this->_FormCovalentBond (atom1, atom2);
	}

	// bonding between N and H
	if ( (atom1->Name().find("N") != string::npos && atom2->Name().find("H") != string::npos)
		|| (atom2->Name().find("N") != string::npos && atom1->Name().find("H") != string::npos) ) {
	
		if (distance < NO_BOND_LENGTH) 
			this->_FormCovalentBond (atom1, atom2);
	}

return;
}

// updates all the fields of the matrix
// 		Top triangle = distance between atoms
// 		Bottom triangle = covalent bonds
void ConnectivityMatrix::UpdateMatrix () {
	
	// first clear out the matrix and resize it appropriately
	_matrix.resize(_atoms.size());

	RUN (_atoms) {
		_matrix[i] = std::vector<double> (_atoms.size(), 0.0);		// zero-out all the elements
		_atoms[i]->ClearHBonds ();								// and let all the atoms also know that we're starting over
	}

	// now run through and calculate the distance between each pair of atoms in the top half of the matrix, and skip the diagonal elements
printf ("starting\n");
	for (unsigned int i=0; i < _atoms.size() - 1; i++) {
		if (_atoms[i]->Residue() != "h2o" && _atoms[i]->Residue().find("no3") == string::npos) continue;
		// let's only worry about water molecules... for now
		for (unsigned int j = i+1; j < _atoms.size(); j++) {
			if (_atoms[j]->Residue() != "h2o" && _atoms[j]->Residue().find("no3") == string::npos) continue;
			
			_matrix[i][j] = *_atoms[i] - *_atoms[j];
			// and form bonds if they exist (covalent/hbonds)
			this->_FindBonds (_atoms[i], _atoms[j]);
		}
	}
printf ("done\n");

return;
}

/* Given a water molecule object, let's use the connectivity matrix to determine its coordination number.
 * What this means is counting up all the H-bond interactions on the hydrogens and the oxygens, and then reporting a number value that corresponds to a given coordination type. */
coordination ConnectivityMatrix::FindWaterCoordination (const Water& water) const {
	
	coordination coord;		// This will be our output
	int donor = 0, acceptor = 0;
	
	for (int atom = 0; atom < water.size(); atom++) {
		// check for donor bonds
		if (water[atom]->Name().find("H") != string::npos) {
			donor += this->HBonds(water[atom]);
		}

		// and for acceptor bonds
		if (water[atom]->Name().find("O") != string::npos) {
			acceptor += this->HBonds(water[atom]);
		}
			
	}
	// copied for reference
	// enum coordination {UNBOUND, O, OH, OHH, OO, OOH, OOHH, H, HH, OOOHH, OOHHH, OVER};

	// now we can run through and set the actual coordination value
	if (donor == 0 && acceptor == 0) coord = UNBOUND;
	else if (donor == 0 && acceptor == 1) coord = O;
	else if (donor == 1 && acceptor == 1) coord = OH;
	else if (donor == 2 && acceptor == 1) coord = OHH;
	else if (donor == 0 && acceptor == 2) coord = OO;
	else if (donor == 1 && acceptor == 2) coord = OOH;
	else if (donor == 2 && acceptor == 2) coord = OOHH;
	else if (donor == 1 && acceptor == 0) coord = H;
	else if (donor == 2 && acceptor == 0) coord = HH;
	else if (donor == 1 && acceptor == 3) coord = OOOH;
	else if (donor == 2 && acceptor == 3) coord = OOOHH;
	else if (donor == 3 && acceptor == 2) coord = OOHHH;
	else if (donor + acceptor > 5) coord = OVER;
	else coord = OTHER;

return coord;
}

std::vector<Atom *> ConnectivityMatrix::ClosestAtoms (const int input, const string atomname, const int number) const {
	
	std::vector< std::vector<double> > distances;	// [distance][test atom ID]
	distances.clear();

	RUN (_atoms) {
		// Skip over the input ID - we're not interested
		if (i == input) continue;
		// also skip over anything that doesn't have the right name
		if (_atoms[i]->Name() != atomname) continue;

		// now we have to store both the test-atom ID and the distance between the test and the input
		std::vector<double> distance;
		distance.push_back (this->Distance(_atoms[input], _atoms[i]));
		distance.push_back (double(i));

		distances.push_back (distance);
	}

	sort (distances.begin(), distances.end());
		
	std::vector<Atom *> output;

	for (int atom = 0; atom < number; atom++) {
		output.push_back (_atoms[int(distances[atom][1])]);
	}

return (output);
}
