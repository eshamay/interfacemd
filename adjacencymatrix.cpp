#include "adjacencymatrix.h"

AdjacencyMatrix::AdjacencyMatrix () :
	_size(0)
{

return;
}

AdjacencyMatrix::AdjacencyMatrix (const Atom_ptr_vec& atoms) :
	_size(atoms.size()) {

	UpdateMatrix (atoms);

return;
}

AdjacencyMatrix::~AdjacencyMatrix () {

	//DeleteMatrix();

return;
}

void AdjacencyMatrix::UpdateMatrix (const Atom_ptr_vec& atoms) {

	// Clear out the old matrix
	//if (_size)
		//DeleteMatrix ();

	// Form the atom list
	_size = atoms.size();
	_atoms.resize (_size, (Atom *)NULL);
	for (unsigned int i = 0; i < _size; i++) {
		_atoms[i] = atoms[i];
	}

	// Build the matrix from the atom set given
	BuildMatrix ();

	// now run through all atom combos and get their bondlengths
	Atom *ai, *aj;

	// this little for-loop bit runs through all atom-pair combinations once and find the bond-types between them
	for (unsigned int i = 0; i < _size-1; i++) {
		for (unsigned int j = i + 1; j < _size; j++) {

			ai = _atoms[i];
			aj = _atoms[j];
			std::string ai_name = ai->Name();
			std::string aj_name = aj->Name();

			// Don't connect oxygens to oxygens, and hydrogen to hydrogen...
			if (ai_name == aj_name) continue;

			// calculate the distance between the two atoms
			double bondlength = ai->Position().MinDistance(aj->Position(), Atom::_size);

			// now set the actual bondtype by checking distance criteria
			bondtype btype = unbonded;

			// first look at bonds between O and H
			if ( (ai_name.find("O") != std::string::npos || aj_name.find("O") != std::string::npos)
				&&
				 (ai_name.find("H") != std::string::npos || aj_name.find("H") != std::string::npos)
			   )
			{

				// one type of bond is the O-H covalent
				if (bondlength <= OHBONDLENGTH) {
					btype = ohbond;
				}

				// Or an H-bond is formed!
				if (bondlength <= HBONDLENGTH && bondlength > OHBONDLENGTH) {
					// additionally, let's check the angle-criteria for an H-bond.
					// This is done by looking at the angle formed from
					#ifdef ANGLE_CRITERIA
					Atom *o1, *h, *o2;	// o1 is covalently bound to h, and o2 is h-bound to h
					Water * wat;

					if (ai->Name().find("O") != std::string::npos) {		// ai is the O, and aj is the H
						o2 = ai;
						h = aj;
					}
					else if (aj->Name().find("O") != std::string::npos) {
						h = ai;
						o2 = aj;
					}
					wat = static_cast<Water *>(h->ParentMolecule());
					o1 = wat->GetAtom("O");

					VecR oh1 = o1->Position().MinVector(h->Position(), Atom::Size());	// the covalent bond
					VecR oh2 = h->Position().MinVector(o2->Position(), Atom::Size());	// the H-bond

					double angle = acos(oh1 < oh2);
					//printf ("% 10.3f\n", angle * 180.0/M_PI);

					if (angle < HBONDANGLE)
					#endif

					btype = hbond;
				}
			}

			// now connect Os to Ns
			if ( (ai_name.find("O") != std::string::npos || aj_name.find("O") != std::string::npos)
				&&
				 (ai_name.find("N") != std::string::npos || aj_name.find("N") != std::string::npos)
			   )
			{
				if (bondlength <= NOBONDLENGTH) {
					btype = nobond;
				}
			}
			// add in the bond between two atoms
			SetBond (i, j, bondlength, btype);
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

	//_matrix.clear();
	_matrix.resize (_size, Bond_vec(_size, Bond()));

	// only set up the upper-triangle of the matrix (don't include the diagonal)
/*
	for (int i = 0; i < _size - 1; i++) {
		for (int j = i + 1; j < _size; j++) {
			_matrix[i][j].bond = unbonded;
		}
	}
*/

return;
}

// Set all the bonds to unbonded
void AdjacencyMatrix::ClearBonds () {

	Bond * b;
	for (unsigned int i = 0; i < _size - 1; i++) {
		for (unsigned int j = i + 1; j < _size; j++) {
			b = &_matrix[i][j];
			b->bond = unbonded;
		}
	}

return;
}

/*
// deletes all the bonds from memory - clean up
void AdjacencyMatrix::DeleteMatrix () {

	Bond * b;

	for (int i = 0; i < _size - 1; i++) {
		for (int j = i + 1; j < _size; j++) {
			if (b == (Bond *)NULL) {
				cout << "\nFound a null bond pointer in the adjacency matrix during AdjacencyMatrix::DeleteMatrix()\n" << endl;
				cout << "i = " << i << " ; j = " << j << endl;
				exit(1);
			}
			b = _matrix[i][j];
			delete b;
		}
	}

	_matrix.clear();
	_atoms.clear();
	_size = 0;

return;
}
*/

void AdjacencyMatrix::SetBond (int x, int y, const double blength, const bondtype btype) {

	// to work with the upper half of the matrix (upper triangle) - fix the indices
	if (x > y) {
		int tmp = x;
		x = y;
		y = tmp;
	}

	Bond * b = &_matrix[x][y];

	b->bondlength = blength;
	b->SetBondType (btype);

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

Atom_ptr_vec AdjacencyMatrix::BondedAtoms (Atom const * const ap) {

	Atom_ptr_vec atoms;

	int id = ID (ap);

	for (int row = 0; row < id; row++) {
		if (_matrix[row][id].bond == unbonded) continue;

		atoms.push_back (_atoms[row]);
//		if (_matrix[row][id]->bond == hbond)
//			printf ("-- %4.3f -->%s) %d\n|",
//				_matrix[row][id]->bondlength, _atoms[row]->Name().c_str(), _atoms[row]->ID());
	}

	for (unsigned int col = id + 1; col < _size; col++) {
		if (_matrix[id][col].bond == unbonded) continue;

		atoms.push_back (_atoms[col]);
//		if (_matrix[id][col]->bond == hbond)
//			/printf ("-- %4.3f -- >%s) %d\n|",
//				_matrix[id][col]->bondlength, _atoms[col]->Name().c_str(), _atoms[col]->ID());
	}

return (atoms);
}

// check if a bond is a covalent bond
bool AdjacencyMatrix::CovalentBond (bondtype const b) const {

	bool cov = false;

	if (b != hbond && b != unbonded)
		cov = true;

	return (cov);
}

// This should return a list of all the atoms hbonded to a given atom (pointer[s])
Atom_ptr_vec AdjacencyMatrix::BondedAtoms (Atom const * const ap, bondtype const b) {

	Atom_ptr_vec atoms;

	int id = ID (ap);

	for (int row = 0; row < id; row++) {
		Bond * b_tmp = &_matrix[row][id];
		bondtype b_type = b_tmp->bond;

		if ( (b == covalent && this->CovalentBond(b_type))
			|| b_type == b)
			atoms.push_back (_atoms[row]);

//		if (_matrix[row][id]->bond == hbond)
//			printf ("-- %4.3f -->%s) %d\n|",
//				_matrix[row][id]->bondlength, _atoms[row]->Name().c_str(), _atoms[row]->ID());
	}

	for (unsigned int col = id + 1; col < _size; col++) {
		Bond * b_tmp = &_matrix[id][col];
		bondtype b_type = b_tmp->bond;

		if ( (b == covalent && this->CovalentBond(b_type))
			|| b_type == b)
			atoms.push_back (_atoms[col]);
//		if (_matrix[id][col]->bond == hbond)
//			/printf ("-- %4.3f -- >%s) %d\n|",
//				_matrix[id][col]->bondlength, _atoms[col]->Name().c_str(), _atoms[col]->ID());
	}

return (atoms);
}

// returns the bound atoms that are named 'name' and bound by a bondtype 'b'.
Atom_ptr_vec AdjacencyMatrix::BondedAtoms (Atom const * const ap, bondtype const b, const std::string name) {

	Atom_ptr_vec atoms;

	int id = ID (ap);

	for (int row = 0; row < id; row++) {
		Bond * b_tmp = &_matrix[row][id];
		bondtype b_type = b_tmp->bond;

		if ( (b == covalent && this->CovalentBond(b_type))
			|| b_type == b) {
			if (_atoms[row]->Name().find(name) != std::string::npos) {
				atoms.push_back (_atoms[row]);
			}
		}

//		if (_matrix[row][id]->bond == hbond)
//			printf ("-- %4.3f -->%s) %d\n|",
//				_matrix[row][id]->bondlength, _atoms[row]->Name().c_str(), _atoms[row]->ID());
	}

	for (unsigned int col = id + 1; col < _size; col++) {
		Bond * b_tmp = &_matrix[id][col];
		bondtype b_type = b_tmp->bond;

		if ( (b == covalent && this->CovalentBond(b_type))
			|| b_type == b) {
			if (_atoms[col]->Name().find(name) != std::string::npos) {
				atoms.push_back (_atoms[col]);
			}
		}
//		if (_matrix[id][col]->bond == hbond)
//			/printf ("-- %4.3f -- >%s) %d\n|",
//				_matrix[id][col]->bondlength, _atoms[col]->Name().c_str(), _atoms[col]->ID());
	}

return (atoms);
}

Bond_ptr_vec AdjacencyMatrix::Bonds (Atom const * const ap) {

	Bond_ptr_vec vb;

	int id = ID (ap);

	for (int row = 0; row < id; row++) {
		if (_matrix[row][id].bond == unbonded) continue;

		vb.push_back (&_matrix[row][id]);
//		if (_matrix[row][id]->bond == hbond)
//			printf ("-- %4.3f -->%s) %d\n|",
//				_matrix[row][id]->bondlength, _atoms[row]->Name().c_str(), _atoms[row]->ID());
	}

	for (unsigned int col = id + 1; col < _size; col++) {
		if (_matrix[id][col].bond == unbonded) continue;

		vb.push_back (&_matrix[id][col]);
//		if (_matrix[id][col]->bond == hbond)
//			/printf ("-- %4.3f -- >%s) %d\n|",
//				_matrix[id][col]->bondlength, _atoms[col]->Name().c_str(), _atoms[col]->ID());
	}

return (vb);
}

// find the bond between two atoms
Bond * AdjacencyMatrix::GetBond (Atom const * const a1, Atom const * const a2) {

	int id1 = ID (a1);
	int id2 = ID (a2);

	// we search the upper triangle of the matrix...
	int temp = id2;
	if (id1 > id2) {
		id2 = id1;
		id1 = temp;
	}

return (&_matrix[id1][id2]);
}

Bond_ptr_vec AdjacencyMatrix::HBonds (Atom const * const ap) {

	Bond_ptr_vec vb = Bonds(ap);
	Bond_ptr_vec hbonds;

	RUN (vb) {
		if (vb[i]->bond == hbond)
			hbonds.push_back(vb[i]);
	}

return (hbonds);
}

int AdjacencyMatrix::NumHBonds (Water const * const wat) {

	int num = 0;
	RUN (wat->Atoms()) {
		num += NumHBonds(wat->Atoms(i));
	}

return (num);
}

// calculates the water bonding coordination of a given water molecule
coordination AdjacencyMatrix::WaterCoordination (Water const * const wat) {

	int c = 0;
	Atom * ap;
	for (int atom = 0; atom < wat->size(); atom++) {
		ap = wat->Atoms(atom);
		int bonds = NumHBonds (ap);


		if (ap->Name().find("H") != std::string::npos) {
			c += 10 * bonds;
		}
		if (ap->Name().find("O") != std::string::npos) {
			c += bonds;
		}
	}
//	printf ("%d)  ", c);

return (coordination) c;
}
