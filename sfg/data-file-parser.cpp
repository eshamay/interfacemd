#include "data-file-parser.h"

int main () {

	datafile_parsers::PolarizabilityDataFile pdf ("h2o_polarizability.dat");

	pdf.Polarizability (0.93, 0.951, 103.4);

	return 0;
}
