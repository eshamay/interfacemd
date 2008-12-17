#include "boxfiller.h"

using namespace std;

BoxFiller::BoxFiller (string paramfile, string pdbfile, double spacing) : _pdb(PDBFile(pdbfile)) {

	// load the parameter file in order to get the box data for populating
	_params = (FILE *)NULL;
	_params = fopen (paramfile.c_str(), "r");
	if (_params == (FILE *)NULL) {
		printf ("Couldn't load the param input file:: %s\n", paramfile.c_str());
		exit(1);
	}

	_residueNames.clear();
	_residueNum.clear();
	_spacing = spacing;

	// grab parameters from the input file
	this->_InitParams();

	// fill up the box!
	this->_FillBoxRandom();

	// Write out the PDB file
	vector <Molecule *> mols;
	for (int i=0; i < _mols.size(); i++) mols.push_back (&_mols[i]);

	PDBFile::WritePDB (mols);

return;
}

BoxFiller::~BoxFiller() { };

void BoxFiller::_FillBoxRandom () {

	Molecule tempMol;

	// For each residue in the input parameter set
	for (int res = 0; res < _residueNames.size(); res++) {

		// we go through the PDB file and find the template residue
		for (int tMol = 0; tMol < _pdb.size(); tMol++) {
			if (_pdb[tMol]->Name() != _residueNames[res]) continue;
			// set the template molecule and prepare to populate the right number of times this particular residue
			tempMol = *_pdb[tMol];
		}

		int numAdded = 0;	// We count the number of residues we've successfully added to the system
		while (numAdded != _residueNum[res]) {
			bool fits = false;				// and whether the new molecule added fits into the system
			Molecule newMol = tempMol;				// create a new molecule to add into the system
			int attempts = 0;

			while (!fits) {

				if (attempts % 100 == 0) { 
					VecR vRot (2.0 * RANDOM - 1.0, 2.0 * RANDOM - 1.0, 2.0*RANDOM-1.0);
					double rot = 2.0 * M_PI * RANDOM - 1.0;
					// randomly rotate it
					newMol.Rotate(VecR(), vRot, rot);
				}

				// randomly move it to somewhere in the box
				VecR shift (RANDOM*_boxSize[x], RANDOM*_boxSize[y], RANDOM*_boxSize[z]);
				newMol.Shift(shift);

				// now check to see if the newly placed molecule fits in the system;
				fits = true;

				for (int mol = 0; mol < numAdded; mol++) {

					if (mol == numAdded) continue;	// for the case of the first molecule added to the system

					double distance = newMol.MinDistance (_mols[mol]);		// here we calculate the distance between the atoms of molecules
					if (distance < _spacing)  {
						//printf ("%f\n", distance);
						fits = false;
						shift *= -1.0;
						newMol.Shift(shift);
						newMol.Rotate(VecR(), vRot, -rot);
						break;
					}
				}
			}

			_mols.push_back(newMol);
			numAdded++;
			printf ("%d  ", numAdded); fflush (stdout);
		}
	}

return;
}

void BoxFiller::_InitParams () {

	double x, y, z;
	int residueNum;
	char residue[1000];

	// get the box size info
	fscanf (_params, " BOX %lf %lf %lf ", &x, &y, &z);
	_boxSize.Set (x,y,z);
	//_boxSize.Print();

	Atom::Size(_boxSize);
	//Atom::Size().Print();
	// get info of all the residues
	while (!feof(_params)) {
		fscanf (_params, " %s %d ", residue, &residueNum);
		_residueNum.push_back (residueNum);
		_residueNames.push_back (string(residue));
	}

	fclose (_params);

	//for (int i = 0; i < _residueNames.size(); i++) {
	//	printf ("%s\t%d\n", _residueNames[i].c_str(), _residueNum[i]);
	//}

return;
}

void BoxFiller::_FillBoxLattice () {

	Molecule tempMol;

	// we need the total number of residues in the system
	int numResidues = _residueNames.size();

	// now calculate the number of rows,columns,layers - we know the volume, so just divide the volume by the num of residues, and find the per-residue amount needed.
	double volume = _boxSize[x] * _boxSize[y] * _boxSize[z];
	double resVolume = volume / numResidues;
	double resSpace = pow(resVolume, 1.0/3.0);

	double shiftX=0.0;
	double shiftY=0.0;
	double shiftZ=0.0;


	while (_residueNames.size()) {
		
		int res = int(RANDOM * _residueNames.size());
		printf ("%d\n", res);
		// we go through the PDB file and find the template residue
		for (int tMol = 0; tMol < _pdb.size(); tMol++) {
			if (_pdb[tMol]->Name() != _residueNames[res]) continue;
			// set the template molecule and prepare to populate
			tempMol = *_pdb[tMol];
		}
	
		Molecule newMol = tempMol;				// create a new molecule to add into the system

		VecR vRot (2.0 * RANDOM - 1.0, 2.0 * RANDOM - 1.0, 2.0*RANDOM-1.0);
		double rot = 2.0 * M_PI * RANDOM - 1.0;
		// randomly rotate the new molecule
		newMol.Rotate(VecR(), vRot, rot);

		VecR shift (shiftX, shiftY, shiftZ);	// a shift vector for the residue in the lattice
		newMol.Shift (shift);					// move the res to the new lattice point

		shiftX += resSpace;
		if (shiftX > _boxSize[x]) {
			shiftX = 0.0;
			shiftY += resSpace;
			if (shiftY > _boxSize[y]) {
				shiftY = 0.0;
				shiftZ += resSpace;
			}
		}

		_mols.push_back (newMol);
		_residueNames.erase(_residueNames.begin() + res);
	}

return;
}
