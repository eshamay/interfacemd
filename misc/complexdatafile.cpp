#include "complexdatafile.h"

// all we need is a path to a valid data file
ComplexDataFile::ComplexDataFile (const char * pathname) {

	_file.open(pathname, std::ios::in);
	if (!_file.is_open()) {
		std::cout << "***** ComplexDataFile::ComplexDataFile - Error - Couldn't open the data file: " << pathname << ".*****\nNow exiting." << std::endl;
		exit (1);
	}
	this->_ParseFile ();

	_file.close();

return;
}

// parsing the file just means taking in all the column data. Each row = 3 columns = dependent variable, real, imaginary.
void ComplexDataFile::_ParseFile () {

	double data, real, imag;

	while (!_file.eof()) {
		_file >> data >> real >> imag;
		_data.push_back (data);
		_complex.push_back (std::complex<double> (real, imag));
	}

return;
}

// do a linear regression to find the data point for the given dependent variable
// The data is organized as 3 columns: Dependent_variable, real data, imaginary data.
std::complex<double> ComplexDataFile::DataPoint (const double data) {

	std::complex<double> value;

	double data_low, data_high;
	int datum_low, datum_high;

	for (unsigned int datum = 0; datum < _data.size(); datum++) {

		// here we find the higher and lower data values
		if (_data[datum] > data) {
			datum_high = datum;
			datum_low = datum - 1;
			break;
		}
	}

	// rise over run.
	double real_slope = (_complex[datum_high].real() - _complex[datum_low].real()) / (_data[datum_high] - _data[datum_low]);
	double imag_slope = (_complex[datum_high].imag() - _complex[datum_low].imag()) / (_data[datum_high] - _data[datum_low]);
	double real_offset = _complex[datum_low].real() - _data[datum_low] * real_slope;
	double imag_offset = _complex[datum_low].imag() - _data[datum_low] * imag_slope;;

	value = std::complex<double> (real_slope * data + real_offset, imag_slope * data + imag_offset);

return value;
}

int main (int argc, char **argv) {

	ComplexDataFile data (argv[1]);

return 0;
}
