//#include "../analysis.h"
//#include "../pdbfile.h"
#include <iostream>
#include <string>

int main () {

	  /*

  PDBFile pdb ("SiO2-slab.pdb");

  Atom * h;
  RUN (pdb.Molecules())
  {
	Molecule * mol = pdb[i];

	vector<Atom *> atoms = mol->Atoms();
	RUN2 (atoms)
	{
	  Atom * atom = atoms[j];

	  VecR pos = atom->Position();

	  if (pos[y] > 13.9)
	  {
		h = new Atom("H", VecR(atom->Position() + VecR(0.0, 1.0, 0.0)));
		mol->AddAtom (h);
	  }
	  else if (pos[y] < 0.1)
	  {
		h = new Atom("H", VecR(atom->Position() - VecR(0.0, 1.0, 0.0)));
		mol->AddAtom (h);
	  }
	}
  }

  PDBFile::WritePDB (pdb.Molecules());
  */
  return 0;
}

/*
#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include "../utility.h"
//#include "../analysis.h"

using namespace std;


int main () {

  // data input from text file
  ifstream input ("mdcrd");
  vector<double> vd;

  istream_iterator<double> ii (input);
  copy (ii, istream_iterator<double>(), back_inserter(vd));
  input.close();

  // data output
  ofstream output ("temp.bin", ios::out | ios::binary);
  output << setprecision(9);
  std::copy (vd.begin(), vd.end(), oi_t<double>(output));
  output.close();

  input.open ("temp.bin", ios::binary);
  ii_t<double> ij(input);
  vector<double> ve;
  copy (ij, ii_t<double>(), back_inserter(ve));
  input.close();

  RUN (ve)
	printf ("%12.7f\n", ve[i]);

  return 0;
}
*/
/*
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
*/

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

/*
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

typedef vector<Atom *>::iterator atom_it;
atom_it atom;
atom_it end = xyz.Atoms().end();

// now rename each atom with a unique identifier
Molecule * mol;
for (atom = xyz.Atoms().begin(); atom != end; atom++)
{
std::string atom_name = (*atom)->Name();
mol = new Molecule();

if (atom_name == "Si") {
mol->Name("qzs");
}
if (atom_name == "O") {
mol->Name("qzo");
}

mol->AddAtom(*atom);
mols.push_back(mol);

}

PDBFile::WritePDB (mols);

return 0;
}

*/

