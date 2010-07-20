#include <vector>
#include <complex>
#include <iostream>
#include <fstream>

class ComplexDataFile {

private:

	std::ifstream				_file;
	std::vector<double>			_data;
	std::vector<std::complex<double> > 	_complex;

	void _ParseFile ();

public:

	ComplexDataFile (const char * path);

	std::complex<double> DataPoint (const double conc);

};
