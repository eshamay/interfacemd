#include "pdbfile.h"

PDBFile::PDBFile (string path) {
  _loaded = 0;

  _path = path;

  _file = (FILE *)NULL;
  _file = fopen (path.c_str(), "r");
  if (_file != (FILE *)NULL)
  {
	_loaded = 1;
  }
  else
  {
	printf ("Couldn't load the PDB coordinate file:: %s\n", path.c_str());
	exit(1);
  }

  // include here a check for the TER
  
  //this->_FindLastStep();

  this->LoadFirst ();	// load the first frame of the file

}

PDBFile::PDBFile () {
}

PDBFile::PDBFile (vector<Molecule *>& mols) {
  _mols = mols;
}

PDBFile::~PDBFile () {
  fclose(_file);
}

void PDBFile::_ParseMolecules()
{
  char line[1000];
  char word[10];

  // first process the frame's header and get to the 'ATOM' entries. We need an initial 'word' to see where we are
  fgets (line, 1000, _file);
  sscanf (line, " %s", word);

  Molecule *pmol;
  pmol = new Molecule;
  //Start parsing atoms until the 'END' of the frame
  while (strcmp(word, "END")) {

	// if we hit an ATOM entry, then parse it into the current working molecule
	if (!strcmp(word, "ATOM")) {
	  // add an atom to the current molecule
	  pmol->AddAtom(_ParseAtom(line));
	  _numAtoms++;
	}

	// The TER keyword separates molecules... so if encountered we need a new one
	if (!strcmp(word, "TER")) {
	  pmol->Name ( (*pmol)[0]->Residue());
	  //printf ("%s\n", pmol->Name().c_str());
	  _mols.push_back(pmol);
	  pmol = new Molecule;
	  _numMols++;
	}

	// then get the next line for processing
	fgets (line, 1000, _file);
	sscanf (line, " %s", word);
  }

  return;
}

void PDBFile::LoadNext () {
  _mols.clear();		// first empty out the list of molecules in the system
  _numAtoms = 0;
  _numMols = 0;

  this->_ParseMolecules();

  _currentstep++;
}

// Take the next line of the file, and parse the
Atom *PDBFile::_ParseAtom (const char *atomEntry) {

  double x, y, z;
  char name[10], residue[10];

  sscanf(atomEntry, "ATOM %*d %s %s %*d %lf %lf %lf", name, residue, &x, &y, &z);
  //printf ("%s\t%s\t%f\t%f\t%f\n", name, residue, x, y, z);

  VecR position (x, y, z);
  string atomName (name);

  Atom *pAtom = new Atom (atomName, position);
  pAtom->Residue (residue);

  return (pAtom);
}

// grab the number of the last time step in the file
int PDBFile::_FindLastStep () {

  rewind(_file);
  _laststep = 0;	// reset the counter

  // PDB files separate frames/timesteps with the END keyword. So running through the file and search for lines that say 'END' will
  // give the number of frames in the file

  char line[1000];
  char word[10];

  while (!feof(_file)) {
	fgets (line, 1000, _file);
	sscanf (line, " %s", word);
	if (!strcmp(word, "END")) {
	  _laststep++;
	}
  }

  rewind(_file);		// rewind the file back to the start to be nice
  _currentstep = 0;

  return (_laststep);
}

void PDBFile::LoadFirst() {
  rewind (_file);
  _currentstep = 0;

  this->LoadNext();

  _atoms.clear();
  RUN (_mols) {
	for (int atom = 0; atom < _mols[i]->size(); atom++) {
	  _atoms.push_back (_mols[i]->Atoms(atom));
	}
  }
}

void PDBFile::Seek (int step) {
  rewind (_file);
  _currentstep = 0;

  while (_currentstep != step) {
	this->LoadNext();
  }
}

void PDBFile::LoadLast () {
  rewind (_file);
  _currentstep = 0;

  while (_currentstep != _laststep) {
	this->LoadNext();
  }
}

void PDBFile::WritePDB (vector<Molecule *>& sys) {

  int atomCount = 1;
  int molCount = 1;

  // go through each molecule
  for (int i = 0; i < (int)sys.size(); i++) {
	Molecule * mol = sys[i];
	// and print out each atom
	for (int atom = 0; atom < mol->size(); atom++) {

	  Atom *tatom = mol->Atoms(atom);
	  printf ("ATOM  %5d  %-4s%3s %5d    % 8.3f% 8.3f% 8.3f\n",
		  atomCount++, tatom->Name().c_str(), tatom->Residue().c_str(), molCount,
		  tatom->X(), tatom->Y(), tatom->Z());
	}

	printf ("TER\n");
	molCount++;
  }

  printf ("END");

  return;
}
