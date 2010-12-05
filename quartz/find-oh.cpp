#include "../pdbfile.h"

/* Find all the oxygens to be connected to the surface hydrogens */

bool TestNamePair (std::string name1, std::string name2);

int main () {

  PDBFile pdb ("temp.pdb");

  Molecule * mol = pdb.Molecules(0);
  mol->Name("qtz");
  std::vector<Atom *> atoms = mol->Atoms();

  /* cycle through each atom in order of the pdb file */
  RUN (atoms)
  {
	Atom * atom = atoms[i];
	atom->Residue("qtz");
	VecR pos = atom->Position();

	/* find the surface Oxygens */
	if (atom->Name().find("O") != string::npos && pos[z] > -7.9)
	{
	  /* add H's to the surface */
	  Atom * h = new Atom ("H", pos + VecR (0.0, 0.0, 1.1));
	  mol->AddAtom(h);
	}
  }
  std::vector<Molecule *> mols;
  mols.push_back(mol);
  PDBFile::WritePDB(mols);
}

bool TestNamePair (std::string name1, std::string name2)
{
  return (name1.c_str()[0] == 'O' && name2.c_str()[0] == 'H') 
	|| 
	(name2.c_str()[0] == 'O' && name1.c_str()[0] == 'H');
}
