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
	void operator() (T t);

	// Output data from the histograms to a file
	void Output (FILE * output, const int timestep);

  private:

  	static vector<string> atom_types;
	// The set of histograms for density binning of atoms - one histogram for each atom-type in the system
	static map<string, Histogram_1D> _histograms;		

	// Create a new histogram in the map
	void AddNewHistogram (T t) {
	  _histograms[t->Name()] = Histogram_1D (Analyzer::posbins, 0);
	  atom_types.push_back(t->Name());
	}

	// Check if the map already contains a histogram for the incoming ... thing
	bool MapKeyExists (T t) {
	  return _histograms.end() != _histograms.find(t->Name());
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
	DensityBinner<Atom *> binner;

	void Setup () { 

	  LoadAll();
	  std::vector<Atom *> v_si;
	  /* grab all the Si atoms */
	  RUN (int_atoms)
	  {
		Atom * atom = int_atoms[i];
		if (atom->Name() == "SI")
		  v_si.push_back(atom);
	  }

	  /* load only the waters */
	  FindWaters();

	  /* now load all the Si atoms into the ones to be processed */
	  copy(v_si.begin(), v_si.end(), back_inserter(int_atoms));

	}

	void Analysis () {
	  // bin each water's position and angle
	  FOR_EACH(int_atoms, binner);
	  return;
	}

	// nothing at the moment
	void PostAnalysis () { return; }

	void DataOutput (const unsigned int timestep) {
	  if (!(timestep % (output_freq * 10)))
		binner.Output(output, timestep);
	  return;
	}
};

#endif
