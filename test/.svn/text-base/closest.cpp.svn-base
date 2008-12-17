#include "../ambersystem.h"
#include "../pdbfile.h"
#include "../utility.h"

//#define mdcrd
#define pdb

int main (int argc, char **argv) {

#ifdef mdcrd
	AmberSystem sys ("prmtop", "mdcrd", "");
#endif
#ifdef pdb
	PDBFile sys (argv[1]);
	Atom::Size(VecR (30., 30., 30.));
#endif

	Atom::Size().Print();
	double closest = 0.0;

	Atom * atom1, * atom2;
	Atom * atomA, * atomB;

	bool first = true;
	RUN (sys) {
	RUN2 (sys) {
		if (i == j) continue;

		atom1 = sys.Atoms(i);
		atom2 = sys.Atoms(j);

		double distance = *atom1 - *atom2;
		if (distance < 0.8) {
			printf ("% 10.4f\n", distance);
			atom1->Print(); atom2->Print();
		}

		if (first || distance < closest) {
			closest = distance;
			atomA = atom1; atomB = atom2;
			if (first) first = false;
		}

	}}

	printf ("closest = % 10.4f\n", closest);
	atomA->Print(); atomB->Print();

return 0;
}
