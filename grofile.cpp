#include "grofile.h"

namespace gromacs {

	GROFile::GROFile (const std::string gropath) :
		md_files::CoordinateFile(gropath)
	{ 
		_ParseSystem();
		return; 
	}


	GROFile::~GROFile () {
		for (Atom_it it = _atoms.begin(); it != _atoms.end(); it++)
			delete ((*it));
		for (Mol_it it = _mols.begin(); it != _mols.end(); it++) 
			delete ((*it));
		return;
	}

	void GROFile::_ParseSystem () {

		rewind (_file);

		ReadLine();		// first line is the header/title
		_title = _line;

		ReadLine();		// clear out the rest of the line
		char temp[100];
		strncpy (temp, &_line[0], 100);
		int num_atoms = atoi(temp);
		
		_atoms.resize(num_atoms, (AtomPtr)NULL);

		int oldresnum = 0;
		int newresnum = 0;

		MolPtr newmol;

		Atom_it_non_const it = _atoms.begin();
		for (int i = 0; i < num_atoms; i++) {
			//std::cout << "i = " << i << std::endl;
			std::string res_id, atom_id;
			std::string resname, atomname;

			// read a line of the file to get all the parameters
			ReadLine();

			// Now parse out all the fields of the GRO file
			// 1-5	Residue number
			res_id.insert(0, &_line[0], 5);
			boost::trim(res_id);
			newresnum = atoi(res_id.c_str());
			//std::cout << "newresnum " << newresnum << std::endl;
			// 6-12 Residue name
			resname.insert(0, &_line[5], 7);
			boost::trim(resname);
			//std::cout << "resname " << resname << std::endl;
			// 13-15 Atom name
			atomname.insert(0, &_line[12], 3);
			boost::trim(atomname);
			//std::cout << "atomname " << atomname << std::endl;
			// 16-20 Atom number
			atom_id.insert(0, &_line[15], 5);
			boost::trim(atom_id);
			//std::cout << "atom_id " << atom_id << std::endl;


			//fscanf (_file, 
			//"%5d%5s%5s%*5d%*8.3f%*8.3f%*8.3f%*8.4f%*8.4f%*8.4f",	// old one
			//"%5d%5s%5s%*5d%*8f%*8f%*8f%*8f%*8f%*8f",		// new one
			//&newresnum, resname, atomname);

			//ReadLine();

			*it = new Atom (std::string(atomname), VecR());
			(*it)->ID(i+1);

			if (newresnum != oldresnum) {
				//std::cout << "making a newmol: old/new = " << oldresnum << "/" << newresnum << std::endl;
				//newmol = new Molecule();
				newmol = molecule::MoleculeFactory(resname);
				//newmol->Name(resname);
				newmol->MolID(newresnum);
				oldresnum = newresnum;
				_mols.push_back(newmol);
			}

			newmol->AddAtom(*it);
			++it;
		}

		ReadLine();	// clear out the rest of the line

		float bx, by, bz;
		sscanf (_line, " %f %f %f ", &bx, &by, &bz);
		_dimensions.Set(bx, by, bz);

		return;
	}

	void GROFile::Print () const {

		printf ("Title: %s", _title.c_str());
		printf ("%d atoms processed\n", (int)_atoms.size());
		printf ("Box dimensions: ");
		_dimensions.Print();
		for (Mol_it it = _mols.begin(); it != _mols.end(); it++) {
			(*it)->Print();
		}

		return;
	}

} // namespace gromacs

/*
	 int main () {

	 GROFile gro ("sw-box.gro");
	 gro.Print();

	 return 0;
	 }
	 */
