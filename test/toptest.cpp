#include "toptest.h"

int main (int argc, char *argv) {

	TOPFile topfile (PRMTOP);

	printf ("number of mols = %d\n", topfile.NumMols());

	for (int i = 0; i < topfile.NumMols(); i++) {
		printf ("%s\n", topfile.MolNames()[i].c_str());
		printf ("%d\n", topfile.MolSizes()[i]);
	}

return 0;
}
