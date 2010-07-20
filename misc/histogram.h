#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_

#include "utility.h"
#include <vector>
#include <iostream>

using namespace std;

template <typename ValueType>
class Histogram {

private:
  
	ValueType _min;
	ValueType _max;
	double _numBins;
	ValueType _binSize;

	ValueType * _data;

public:

	Histogram (ValueType const min, ValueType const max, ValueType const binsize);
	~Histogram ();

	ValueType Min () const { return _min; }
	ValueType Max () const { return _max; }
	int NumBins () const { return _numBins; }
	ValueType BinSize () const { return _binSize; }

	int Bin (ValueType const value) const {		// Returns the bin that the value would go into
		return int( (value - _min)/_binSize );
	}

	int Add (ValueType const value) 			// Increment the histogram bin for the given value
	{
		int bin = Bin(value);
		_data[bin]++;
		return bin;
	}
	int Remove (ValueType const value)			// Decrement the associated bin for the value given
	{
		int bin = Bin(value);
		_data[bin]--;
		return bin;
	}

};

template <typename ValueType>
Histogram::Histogram (ValueType const min, ValueType const max, ValueType const binsize)	
	: _min (min), _max (max), _binSize (binsize), _numBins ( int((max - min)/binsize) + 1)
{
	_data = new ValueType [_numBins];
return;
}
	
<typename ValueType>
Histogram::~Histogram () {
	delete [] _data;
return;
}

<typename ValueType>
void Histogram::Histogram Clear () {
	for (int i = 0; i < _numBins; i++) {
		_data[i] = ValueType(0);
	}
return;
}

#endif
