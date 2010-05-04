#ifndef __DENSITY_HISTOGRAM_H
#define __DENSITY_HISTOGRAM_H

#include "../utility.h"
#include "../analysis.h"
#include <libconfig.h++>

typedef vector<int> Histogram_1D;
typedef map<string, Histogram_1D>::iterator Histogram_it;


/********************************************************************/
// A binning functor for accumulating a 2D histogram of angles and positions of molecules
// The functor can operate on a molecule so long as the molecule has a method for returning its molecular axis (vecr)
template <class T, class U>
class DensityBinner {
  public:

	// The functor's operation for binning of atomic positions
	void operator() (T t);

	// Output data from the histograms to a file
	void Output (FILE * output, const int timestep);

  private:

  	static vector<string> atom_types;
	// The set of histograms for density binning of atoms - one histogram for each atom-type in the system
	static map<string, Histogram_1D> _histograms;		

	// Create a new histogram in the map
	void AddNewHistogram (T t) {
	  _histograms[t->Name()] = Histogram_1D (Analyzer<U>::posbins, 0);
	  atom_types.push_back(t->Name());
	}

	// Check if the map already contains a histogram for the incoming ... thing
	bool MapKeyExists (T t) {
	  return _histograms.end() != _histograms.find(t->Name());
	}
};

template <class T, class U> vector<string> DensityBinner<T,U>::atom_types = vector<string>();
template <class T, class U> map<string, Histogram_1D> DensityBinner<T,U>::_histograms;


/********************************************************************/
/********************************************************************/
/********************************************************************/


template<class T>
class DensityAnalyzer : public Analyzer<T> {

  public:
	DensityAnalyzer (WaterSystemParams& wsp) :
	  Analyzer<T>(wsp) { return; }

	// the utility functor for getting all the data accumulated
	DensityBinner<Atom *, T> binner;

	void Setup () { 
	  /* Load all the atoms in the system */
	  this->LoadAll();
	}

	void Analysis () {
	  // bin each atom's position in the system
	  std::for_each(this->int_atoms.begin(), this->int_atoms.end(), binner);
	  return;
	}

	// nothing at the moment
	void PostAnalysis () { return; }

	void DataOutput (const unsigned int timestep) {
	  if (!(timestep % (this->output_freq * 10)))
		binner.Output(this->output, timestep);
	  return;
	}
};

#endif
