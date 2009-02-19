#include "adjacencymatrix.h"

AdjacencyMatrix::AdjacencyMatrix () : _size(0) {

return;
}

AdjacencyMatrix::AdjacencyMatrix (const Atom_ptr_vec& atoms) : 
	_size(atoms.size()) {

	UpdateMatrix (atoms);

return;
}

AdjacencyMatrix::~AdjacencyMatrix () {

	DeleteMatrix();

return;
}

void AdjacencyMatrix::UpdateMatrix (const Atom_ptr_vec& atoms) {
	

	// Clear out the old matrix
	if (_size)
		DeleteMatrix ();

	// Form the atom list
	_size = atoms.size();
	_atoms.resize (_size, (Atom *)NULL);
	for (int i = 0; i < _size; i++) {
		_atoms[i] = atoms[i];
	}

	// Build the matrix from the atom set given
	BuildMatrix ();

	// now run through all atom combos and get their bondlengths
	double bondlength = 0.0;
	Atom *ai, *aj;
	
	// this little for-loop bit runs through all atom-pair combinations once and find the bond-types between them
	for (unsigned int i = 0; i < _size-1; i++) {
		for (unsigned int j = i + 1; j < _size; j++) {
			
			ai = _atoms[i];
			aj = _atoms[j];

			// Don't connect oxygens to oxygens, and hydrogen to hydrogen...
			if (
				(ai->Name().find("H") != string::npos and aj->Name().find("H") != string::npos)
				or
				(ai->Name().find("O") != string::npos and aj->Name().find("O") != string::npos)
			)
				continue;
			
			// calculate the distance between the two atoms
			//double distance = *_atoms[i] - *_atoms[j];
			//double distance = _atoms[i]->MinDistance(_atoms[j]);
			double distance = ai->Position().MinDistance(aj->Position(), Atom::_size);

			// check that the distance is at least an H-bond
			if (distance > HBONDLENGTH) continue;

			// add in the bond between two atoms
			SetBond (i, j, distance);
		}
	}

return;
}

// Allocates memory and makes all the Bonds available.
// Make sure the _size is already set before calling this
void AdjacencyMatrix::BuildMatrix () {

	if (!_size) {
		cout << "\nAdjacencyMatrix::BuildMatrix () - trying to build a matrix of size 0" << endl;
		exit(1);
	}

	_matrix.resize (_size, Bond_ptr_vec(_size, (Bond *)NULL));
	
	// only set up the upper-triangle of the matrix (don't include the diagonal)
	for (int i = 0; i < _size - 1; i++) {
		for (int j = i + 1; j < _size; j++) {
			_matrix[i][j] = new Bond;
			_matrix[i][j]->bond = unbonded;
		}
	}

return;
}

// Set all the bonds to unbonded
void AdjacencyMatrix::ClearBonds () {

	Bond * b;
	for (int i = 0; i < _size - 1; i++) {
		for (int j = i + 1; j < _size; j++) {
			b = _matrix[i][j];
			b->bond = unbonded;
		}
	}
	
return;
}

// deletes all the bonds from memory - clean up
void AdjacencyMatrix::DeleteMatrix () {

	Bond * b;
	for (int i = 0; i < _size - 1; i++) {
		for (int j = i + 1; j < _size; j++) {
			b = _matrix[i][j];
			delete b;
		}
	}

	_matrix.clear();
	_atoms.clear();
	_size = 0;

return;
}

void AdjacencyMatrix::SetBond (int x, int y, const double length) const {

	// to work with the upper half of the matrix (upper triangle) - fix the indices
	if (x > y) {
		int tmp = x;
		x = y;
		y = tmp;
	}

	Bond * b = _matrix[x][y];

	b->bondlength = length;
	b->SetBondType ();

return;
}

// returns the atom's location in the matrix
int AdjacencyMatrix::ID (Atom const * const ap) const {
	
	int id = -1;

	for (unsigned int i = 0; i < _size; i++) {
		if (_atoms[i] == ap)
			id = i;
	}
		
	if (id == -1) {
		cout << "\nAdjacencyMatrix::ID (Atom * ap)\nCouldn't find the atom in the matrix" << endl;
		exit (1);
	}

return (id);
}

Bond_ptr_vec AdjacencyMatrix::Bonds (Atom const * const ap) const {

	Bond_ptr_vec vb;

	int id = ID (ap);

	for (int row = 0; row < id; row++) {
		if (_matrix[row][id]->bond == unbonded) continue;

		vb.push_back (_matrix[row][id]);
//		if (_matrix[row][id]->bond == hbond) 
//			printf ("-- %4.3f -->%s) %d\n|", 
//				_matrix[row][id]->bondlength, _atoms[row]->Name().c_str(), _atoms[row]->ID());
	}

	for (int col = id + 1; col < _size; col++) {
		if (_matrix[id][col]->bond == unbonded) continue;

		vb.push_back (_matrix[id][col]);
//		if (_matrix[id][col]->bond == hbond) 
//			/printf ("-- %4.3f -- >%s) %d\n|", 
//				_matrix[id][col]->bondlength, _atoms[col]->Name().c_str(), _atoms[col]->ID());
	}

return (vb);
}

int AdjacencyMatrix::NumBonds (Atom const * const ap) const {

return (Bonds(ap).size());
}

int AdjacencyMatrix::NumHBonds (Atom const * const ap) const {

//********
//printf ("%s) %d\n|", ap->Name().c_str(), ap->ID());

	Bond_ptr_vec vb = Bonds (ap);

	int num = 0;
	for (int i = 0; i < vb.size(); i++) {
		if (vb[i]->bond == hbond)
			num++;
	}

return (num);
}

// calculates the water bonding coordination of a given water molecule
coordination AdjacencyMatrix::WaterCoordination (Water const * const wat) const {
	
	int c = 0;
	Atom * ap;
	for (int atom = 0; atom < wat->size(); atom++) {
		ap = wat->Atoms(atom);
		int bonds = NumHBonds (ap);


		if (ap->Name().find("H") != string::npos) {
			c += 10 * bonds;
		}
		if (ap->Name().find("O") != string::npos) {
			c += bonds;
		}
	}
//	printf ("%d)  ", c);

return (coordination) c;
}
