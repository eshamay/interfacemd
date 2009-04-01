#include "../ambersystem.h"
#include "../pdbfile.h"
#include "../utility.h"

using namespace std;

int main (int argc, char **argv) {

	AmberSystem sys ("prmtop", "mdcrd", "");

	PDBFile::WritePDB (sys.Molecules());

return 0;
}
