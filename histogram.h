#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_

#include "utility.h"
#include <vector>
#include <iostream>

using namespace std;

class Histogram {

private:
	
	double _min;
	double _max;
	double _numBins;
	double _binSize;

	vector<int> _histogram;

public:
	
	Histogram (double const min, double const max, int numbins);
	Histogram (double const min, double const max, double binsize);

	double Min () { return _min; }
	double Max () { return _max; }
	double NumBins () { return _numBins; }
	double BinSize () { return _binSize; }

	int		operator[] (int bin) { return _histogram[bin]; }
	void	operator[] (double value) { _histogram[this->Bin(value)]++; }
	int 	operator++ (int bin) { _histogram[bin]++; }
	int 	operator-- (int bin) { _histogram[bin]--; }

	int Bin (int const value) const {		// calculates the bin for the given value
		return int( (value - _min)/_binSize );
	}

	void Print ();
};

#endif
