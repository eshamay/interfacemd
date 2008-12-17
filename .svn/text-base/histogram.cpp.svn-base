#include "histogram.h"

Histogram::Histogram (double const min, double const max, int numbins) {

	_min = min;
	_max = max;
	_numBins = numbins;
	_binSize = (_max - _min) / _numBins;

	_histogram.clear();
	_histogram.resize (_numBins, 0);
}

Histogram::Histogram (double const min, double const max, double binsize) {
	
	_min = min;
	_max = max;
	_binSize = binsize;
	_numBins = (_max - _min) / _binSize;

	_histogram.clear();
	_histogram.resize (_numBins, 0);
}

void Histogram::Print () {
	
	for (int i = 0; i < _numBins; i++)
		
		printf ("% 8.4f\t% d\n", _min + double(i)*_binSize, _histogram[i]);

return;
}
