#ifndef __DENSITY_HISTOGRAM_H
#define __DENSITY_HISTOGRAM_H

#include "../utility.h"
#include "../analysis.h"

typedef vector<int> Histogram_1D;
typedef map<string, Histogram_1D>::iterator Histogram_it;


/********************************************************************/
// A binning functor for accumulating a 2D histogram of angles and positions of molecules
// The functor can operate on a molecule so long as the molecule has a method for returning its molecular axis (vecr)
template <class T>
class DensityBinner {
  public:

	// The functor's operation for binning of atomic positions
	void operator() (T * t);

	// Output data from the histograms to a file
	void Output (FILE * output);

  private:

  	static vector<string> atom_types;
	// The set of histograms for density binning of atoms - one histogram for each atom-type in the system
	static map<string, Histogram_1D> _histograms;		

	// Create a new histogram in the map
	void AddNewHistogram (T * t) {
	  _histograms[t->Name()] = Histogram_1D (Analyzer::posbins, 0);
	  atom_types.push_back(t->Name());
	}

	// Check if the map already contains a histogram for the incoming ... thing
	bool MapKeyExists (T * t) {
	  Histogram_it it = _histograms.find(t->Name());
	  return (it != _histograms.end());
	}
};

template <class T> vector<string> DensityBinner<T>::atom_types = vector<string>();
template <class T> map<string, Histogram_1D> DensityBinner<T>::_histograms;


/********************************************************************/
/********************************************************************/
/********************************************************************/


class DensityAnalyzer : public Analyzer {

  public:
	DensityAnalyzer (WaterSystemParams& wsp) :
	  Analyzer(wsp) { return; }

	// the utility functor for getting all the data accumulated
	DensityBinner<Atom> binner;

	void Setup () { LoadAll(); }

	void Analysis () {
	  // bin each water's position and angle
	  FOR_EACH(int_atoms, this->binner);
	  return;
	}

	// nothing at the moment
	void PostAnalysis () { return; }

	void DataOutput (const unsigned int timestep) {
	  if (!(timestep % (output_freq * 10)))
		this->binner.Output(output);
	  return;
	}
};

#endif
