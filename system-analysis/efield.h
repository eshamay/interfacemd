#ifndef __EFIELD_H
#define __EFIELD_H

#include "../utility.h"
#include "../analysis.h"
#include <libconfig.h++>

typedef vector<double> Histogram_1D;
typedef Histogram_1D::const_iterator histo_1d_it;

/********************************************************************/
// A binning functor for accumulating a 2D histogram of angles and positions of molecules
// The functor can operate on a molecule so long as the molecule has a method for returning its molecular axis (vecr)
template <class U>
class EFieldBinner {
	public:

		EFieldBinner () {
			_histogram.resize(Analyzer<U>::posbins, 0.0);
			return;
		}

		// The functor's operation for binning of atomic positions
		void operator() (Atom * t);

		// Output data from the histograms to a file
		void Output (FILE * output, const int timestep);

	private:

		static Histogram_1D _histogram;
		static double cutoff;
		double AtomCharge (Atom * t);
		double efield, distance, position, charge;
		int bin;
};

template <class U> Histogram_1D EFieldBinner<U>::_histogram;
template <class U> double EFieldBinner<U>::cutoff = 10.0;

/********************************************************************/
/********************************************************************/
/********************************************************************/

template<class T>
class EFieldAnalyzer : public Analyzer<T> {

	public:
		EFieldAnalyzer (WaterSystemParams& wsp) :
			Analyzer<T>(wsp) { 
				WaterSystem<T>::posmin = wsp.config_file->lookup("analysis.efield.position-cutoff-low");
				WaterSystem<T>::posmax = wsp.config_file->lookup("analysis.efield.position-cutoff-high");
				printf ("\n\tLimiting the E-field analysis to atoms within the region : % 8.3f - %8.3f\n\n",
						WaterSystem<T>::posmin, WaterSystem<T>::posmax);
				return; 
			}

		// the utility functor for getting all the data accumulated
		EFieldBinner<T> binner;

		void Setup () { 
			/* Load all the atoms in the system */
			this->LoadAll();
			std::pair<double,double> extent = std::make_pair(WaterSystem<T>::posmin, WaterSystem<T>::posmax);
			this->KeepAtomsInSlice(this->int_atoms, extent);
		}

		void Analysis () {
			// bin each atom's position in the system
			Atom_ptr_vec& atoms = this->int_atoms;
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
