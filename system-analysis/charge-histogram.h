#ifndef __CHARGE_HISTOGRAM_H
#define __CHARGE_HISTOGRAM_H

#include "../utility.h"
#include "../analysis.h"
#include <libconfig.h++>

typedef vector<double> Histogram_1D;
typedef Histogram_1D::const_iterator histo_1d_it;

/********************************************************************/
// A binning functor for accumulating a 2D histogram of angles and positions of molecules
// The functor can operate on a molecule so long as the molecule has a method for returning its molecular axis (vecr)
template <class U>
class ChargeBinner {
  public:

    ChargeBinner () {
      _histogram.resize(Analyzer<U>::posbins, 0.0);
      return;
    }

    // The functor's operation for binning of atomic positions
    void operator() (Atom * t);

    // Output data from the histograms to a file
    void Output (FILE * output, const int timestep);

  private:

    static Histogram_1D _histogram;
    double AtomCharge (Atom * t);
};

template <class U> Histogram_1D ChargeBinner<U>::_histogram;

/********************************************************************/
/********************************************************************/
/********************************************************************/

template<class T>
class ChargeAnalyzer : public Analyzer<T> {

  public:
    ChargeAnalyzer (WaterSystemParams& wsp) :
      Analyzer<T>(wsp) { return; }

    // the utility functor for getting all the data accumulated
    ChargeBinner<T> binner;

    void Setup () { 
      /* Load all the atoms in the system */
      this->LoadAll();
    }

    void Analysis () {
      // bin each atom's position in the system
      Atom_ptr_vec& atoms = this->sys_atoms;
      std::for_each(atoms.begin(), atoms.end(), binner);

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
