#include "grofile.h"

void GROFile::_ParseSystem () {

  rewind (_file);
  char title[500];
  fgets (title, 500, _file);		// first line is the header/title
  _title = title;
  fscanf (_file, " %d", &_natoms);	// grab the number of atoms
  fgets (title, 500, _file);		// clear out the rest of the line

  _atoms.resize(_natoms, (Atom *)NULL);

  int oldresnum = 1;
  int newresnum = 1;
  char resname[5];
  char atomname[5];
  Molecule * newmol = new Molecule();
  Atom * newatom;

  for (int i = 0; i < _natoms; i++) {
    // read a line of the file to get all the parameters
    fscanf (_file, 
	"%5d%5s%5s%*5d%*8.3f%*8.3f%*8.3f%*8.4f%*8.4f%*8.4f", 
	&newresnum, resname, atomname);
    fgets (title, 500, _file);
    _atoms[i] = new Atom (std::string(atomname), VecR());
    _atoms[i]->ID(i+1);

    if (newresnum != oldresnum) {
      _mols.push_back(newmol);
      newmol = new Molecule();
    }
    newmol->Name(resname);
    newmol->MolID(newresnum);
    oldresnum = newresnum;
    newmol->AddAtom(_atoms[i]);
  }
  _mols.push_back(newmol);
  fgets (title, 500, _file);

  int bx, by, bz;
  fscanf (_file, " %f %f %f ", &bx, &by, &bz);
  _box_size.Set(bx, by, bz);

  return;
}

void GROFile::Print () const {

  printf ("Title: %s", _title.c_str());
  printf ("%d atoms processed\n", _natoms);
  printf ("Box dimensions: ");
  _box_size.Print();
  for (Mol_it it = _mols.begin(); it != _mols.end(); it++) {
    (*it)->Print();
  }

  return;
}

/*
int main () {

  GROFile gro ("sw-box.gro");
  gro.Print();

  return 0;
}
*/
