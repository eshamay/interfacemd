#include "../pdbfile.h"
#include "../mdsystem.h"


int main () {

  PDBFile pdb ("octadecane.pdb");
  //MDSystem::Dimensions(VecR(29.868, 100.0, 27.792));

  Atom_ptr_vec& atoms = pdb.Atoms();
	for (Atom_it it = atoms.begin(); it != atoms.end(); it ++) {
		(*it)->Print();
	}

  printf ("!entry.qtz.unit.connectivity table  int atom1x  int atom2x  int flags\n");
  /* cycle through each atom in order of the pdb file */
  for (int i = 0; i != atoms.size() - 1; i++)
  {
	for (int j = i+1; j != atoms.size(); j++)
	{
	  //double distance = MDSystem::Distance(atoms[i],atoms[j]).Magnitude();
	  double distance = (atoms[i]->Position() - atoms[j]->Position()).Magnitude();
	  if (distance < 1.70)
	  {
		printf (" %d %d 1\n", i+1, j+1);
	  }
	}
  }

  return 0;
}
