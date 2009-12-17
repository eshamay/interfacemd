#include "../pdbfile.h"

int main () {

  PDBFile pdb ("SiO2-slab+H.independent-sio2-atoms.pdb");

  typedef std::vector<Atom *> va;
  typedef std::vector<Molecule *> vm;

  va& atoms = pdb.Atoms();
  vm mols;

  RUN (atoms)
  {
	Atom * atom = atoms[i];

	// create a new residue for each atom
	Molecule * mol = new Molecule();

	// set the residue name 
	std::string name;
	if (atom->Name() == "SI")
	  name = "qsi";
	else if (atom->Name() == "O")
	  name = "qso";
	mol->Name(name);
	mols.push_back(mol);

	// add the atom into the current molecule
	mol->AddAtom(atom);

	/* now add an H to the outside of the SiO2 slab surfaces */
	VecR pos = atom->Position();

	if (pos[y] > 8.9)
	{
	  // rename the O bound to the H
	  mol->Name("qsb");
	  atom->Residue("qsb");
	  // and add an H residue on the slab surface
	  Atom * h = new Atom("H", pos + VecR(0.0, 1.0, 0.0));
	  mol->AddAtom(h);
	}
	else if (pos[y] < 0.1)
	{
	  // rename the O bound to the H
	  mol->Name("qsb");
	  atom->Residue("qsb");
	  // and add an H residue on the slab surface
	  Atom * h = new Atom("H", pos - VecR(0.0, 1.0, 0.0));
	  mol->AddAtom(h);
	}
  }

  PDBFile::WritePDB (mols);
  return 0;
}
