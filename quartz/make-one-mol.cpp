#include "../pdbfile.h"

std::string itoa(int value);

int main () {

  PDBFile pdb ("SiO2-slab.one-res.pdb");

  typedef std::vector<Atom *> va;
  typedef va::iterator vai;

  // the number of Si, H, and O atoms processed - to keep track of the names
  int si = 0;
  int o = 0;
  int h = 0;

  va atoms; // the new container for atoms

  vai va_end = pdb.Atoms().end();
  for (vai atom = pdb.Atoms().begin(); atom != va_end; atom++)
  {
	atoms.push_back(*atom);

	int new_value;
	// change the name of the atom to something more appropriate
	std::string name = (*atom)->Name();
	std::string new_name ("");

	if (name.find("SI") != std::string::npos) {
	  new_name = "S";
	  new_value = ++si;
	} 
	else if (name.find("O") != std::string::npos) {
	  new_name = "O";
	  new_value = ++o;
	}
	new_name += itoa(new_value / 26) + itoa(new_value % 26);
	(*atom)->Name(new_name);




	/* now add an H to the outside of the SiO2 slab surfaces */
	VecR pos = (*atom)->Position();

	if (pos[y] > 8.9)
	{
	  new_value = ++h;
	  new_name = "H";
	  new_name += itoa(new_value / 26) + itoa(new_value % 26);
	  Atom * h = new Atom(new_name, pos + VecR(0.0, 1.0, 0.0));
	  atoms.push_back(h);
	}
	else if (pos[y] < 0.1)
	{
	  new_value = ++h;
	  new_name = "H";
	  new_name += itoa(new_value / 26) + itoa(new_value % 26);
	  Atom * h = new Atom(new_name, pos - VecR(0.0, 1.0, 0.0));
	  atoms.push_back(h);
	}
  }

  Molecule * mol = new Molecule();
  mol->Name("qtz");
  va_end = atoms.end();
  for (vai atom = atoms.begin(); atom != va_end; atom++)
  {
	mol->AddAtom(*atom);
  }
  std::vector<Molecule *> mols;
  mols.push_back(mol);

  PDBFile::WritePDB (mols);
  return 0;
}




/**
 * C++ version std::string style "itoa":
 */

std::string itoa(int value) {
  const char digitMap[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  std::string buf;

  // Guard:
  if (value < 0 || value > 26) {
	// Error: may add more trace/log output here
	std::cout << "std::string itoa(int value) --> error, value must be > 0 and < 26" << std::endl;
	return buf;
  }

  // Translating number to string with base:
  buf = digitMap[value];
  return buf;
}
