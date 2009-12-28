//#include "../watersystem.h"
#include "../xyzsystem.h"
#include "../pdbfile.h"

using namespace std;

int main () {
  XYZFile xyz ("beta-cristobalite-unitcell.xyz");
 
  // create a molecule and add in all the atoms
  Molecule mol;
  mol.name("qtz");
  RUN (xyz)
	mol.AddAtom(xyz[i]);

// now rename each atom with a unique identifier
for (int i = 1; i < 26; i++) 
  	std::cout << itoa(i) << std::endl;

return 0;
}
