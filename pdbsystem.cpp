#include "pdbsystem.h"

namespace md_files {

	void PDBSystem::CountAtoms () {
		_numAtoms = 0;

		rewind(_file);
		char word[6];

		do {
			ReadLine();
			sscanf (_line, "%6s", word);
			if (!strcmp(word, "ATOM"))
				_numAtoms++;
		} while (strcmp(word, "END"));
	}

	void PDBSystem::CreateAtoms () {
		if (!_initialized) {
			for (Atom_it it = begin(); it != end(); it++)
				delete *it;

			_atoms.clear();
			for (int i = 0; i < _numAtoms; i++)
				_atoms.push_back(new Atom);

		}
	}


	void PDBSystem::_ParseAtoms() {
		char keyword[6], atom_id[5], molid[4], X[8], Y[8], Z[8], temp[10];

		for (Atom_it it = begin(); it != end(); it++) {
			// skip all non-ATOM lines
			do {
				ReadLine();
				sscanf (_line, "%6s", keyword);
			} while (strcmp(keyword, "ATOM"));

			std::string atom_name, residue, location;

			// PDB files are fixed-width field entries as described in: http://bmerc-www.bu.edu/needle-doc/latest/atom-format.html

			// 1-6	"ATOM  "
			strncpy (temp, &_line[0], 6);
			// 7-11		atom id number
			strncpy (atom_id, &_line[6], 6);
			// 12 whitespace
			strncpy (temp, &_line[11], 1);
			// 13-16 atom name
			atom_name.insert(0, &_line[12], 4);
			//strncpy (&atom_name[0], &_line[12], 4);
			// 17 location
			//strncpy (&location[0], &_line[16], 1);
			location.insert(0, &_line[16], 1);
			// 18-20 residue name
			//strncpy (&residue[0], &_line[17], 3);
			residue.insert(0, &_line[17], 3);
			// 21 whitespace
			strncpy (temp, &_line[20], 1);
			// 22 chain identifier
			strncpy (temp, &_line[21], 1);
			// 23-26 residue id
			strncpy (molid, &_line[22], 4);
			// 27-30	misc (?) 
			strncpy (temp, &_line[26], 4);
			// 31-38	x-coordinate
			strncpy (X, &_line[30], 8);
			// 39-48	y-coordinate
			strncpy (Y, &_line[38], 8);
			// 47-54	z-coordinate
			strncpy (Z, &_line[46], 8);


			VecR position (atof(X), atof(Y), atof(Z));
			boost::trim(atom_name);
			(*it)->Name(atom_name);
			(*it)->Position(position);
			(*it)->ID(atoi(atom_id));
			(*it)->MolID(atoi(molid));
			(*it)->Residue (residue);
			(*it)->SetAtomProperties();
		}
	}


	void PDBSystem::_ParseMolecules() {
		for (Mol_it it = begin_mols(); it != end_mols(); it++)
			delete *it;
		_mols.clear();

		MolPtr pmol = (MolPtr)NULL;
		int _currentmol = 0;

		for (Atom_it it = begin(); it != end(); it++) {

			// then, if the atom's molid is different than the one we were just on, then readjust the current molecule to something else
			if ((*it)->MolID() != _currentmol) {

				MolPtr newmol;
				newmol = molecule::MoleculeFactory ((*it)->Residue());

				_currentmol = (*it)->MolID();
				newmol->MolID (_currentmol);
				newmol->Name ((*it)->Residue());

				if (pmol != (MolPtr)NULL) {
					_mols.push_back (pmol);
				}

				pmol = newmol;
			}

			pmol->AddAtom(*it);
		}
	}	// Parse Molecules



	void PDBSystem::LoadNext () {
		_ParseAtoms();
		_ParseMolecules();
	}


	//ATOM      1  H1  h2o A   1      27.538  20.316  20.655  1.00  0.00        
	//ATOM    145  N   VAL A  25      32.433  16.336  57.540  1.00 11.92      A1   N
	//123456789a123456789b123456789c123456789d123456789e123456789f
  //         1         2         3         4         5         6         7         8 
	//sscanf(pdbline, "ATOM  %5d %4c%c%3c %c%4d%c%8.3lf%8.3lf%8.3lf", &atomid, name, &location, residue, temp, &molid, temp, &x, &y, &z);
	//sscanf(pdbline, "ATOM  %5d %4c%3c %c%4d%s%8lf%8lf%8lf", &atomid, name, residue, temp, &molid, temp, &X, &Y, &Z);

	void PDBSystem::WritePDB (const Mol_ptr_vec& mols, FILE * file) {

		int atomCount = 1;
		int molCount = 1;

		// go through each molecule
		for (Mol_it it = mols.begin(); it != mols.end(); it++) {

			// and print out each atom
			for (Atom_it jt = (*it)->begin(); jt != (*it)->end(); jt++) {

				fprintf (file, "ATOM  %5d  %-4s%3s %5d    % 8.3f% 8.3f% 8.3f\n",
						atomCount++, (*jt)->Name().c_str(), (*jt)->Residue().c_str(), molCount,
						(*jt)->X(), (*jt)->Y(), (*jt)->Z());
			}

			fprintf (file, "TER\n");
			molCount++;
		}

		fprintf (file, "END");

		return;
	}

}	// namespace md_files
