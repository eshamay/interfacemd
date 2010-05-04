#include <iostream>
#include "../xyzsystem.h"
#include "../pdbfile.h"
#include <map>

using namespace std;

std::string itoa(int value);

int main () {
  XYZFile xyz ("SiO2-slab.xyz");

  vector<Molecule *> mols;
  // create a molecule and add in all the atoms
  Molecule mol;
  mol.Name("qtz");
  RUN (xyz)
	mol.AddAtom(xyz[i]);

  mols.push_back(&mol);

  typedef vector<Atom *>::iterator atom_it;
  atom_it atom;
  atom_it end = xyz.Atoms().end();

  typedef map<std::string, int> name_map;
  name_map names;

  // now rename each atom with a unique identifier
  for (atom = xyz.Atoms().begin(); atom != end; atom++)
  {
	std::string atom_name = (*atom)->Name();
	if (atom_name == "Si")
	  atom_name = "S";
	
	// if the name doesn't already exist in the map
	if (names.find(atom_name) == names.end())
	{
	  // add it in
	  names[atom_name] = 0;
	}

	int value = names[atom_name];
	std::string new_name = atom_name + itoa(value / 26) + itoa(value % 26);
	names[atom_name]++;

	(*atom)->Name(new_name);
  }

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

