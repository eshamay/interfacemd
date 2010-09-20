#include "data-file-parser.h"

int main () {

	datafile_parsers::DipoleDataFile ddf ("h2o_dipole.dat");

	ddf.Value(0.974, 1.532, 105.4).Print();

	return 0;
}
