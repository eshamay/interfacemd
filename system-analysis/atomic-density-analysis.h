#ifndef __DENSITY_HISTOGRAM_H
#define __DENSITY_HISTOGRAM_H

#include "../utility.h"
#include "../analysis.h"
//#include <libconfig.h++>

//typedef vector<int> Histogram_1D;
//typedef map<string, Histogram_1D>::iterator Histogram_it;

template <class T>
class atomic_density_analysis : public AnalysisSet< Analyzer<T> > {
	public:
		typedef Analyzer<T> system_t;

		atomic_density_analysis () : 
			AnalysisSet<system_t> (
					std::string("An analysis of the density of atoms in a system based on atomic position"),
					std::string("3d-atomic-density.dat")
					)
			{ }

		void Setup (system_t& t);
		void Analysis (system_t& t);
		// For each atom type (name) in the system, the histograms in each direction will be output
		void DataOutput (system_t& t);

	protected:
		std::vector<std::string> atom_name_list;
		// Every atom-name will have its own histogram of positions in the system. Each position is held as a vector to the atom site.
		typedef histogram_utilities::Histogram1D<double>	histogram_t;
		typedef std::vector<histogram_t>									histogram_set;
		typedef std::pair<std::string, histogram_set>			histogram_map_elmt;
		typedef std::map<std::string, histogram_set>			histogram_map;
		histogram_map histograms;

		class atomic_position_binner : public std::binary_function<AtomPtr,histogram_map *,void> {
			public:
				void operator() (AtomPtr atom, histogram_map * histos) const {
					histogram_set& hs = histos->operator[] (atom->Name());
					hs[0].operator() (atom->X());
					hs[1].operator() (atom->Y());
					hs[2].operator() (atom->Z());
				}
		} binner;

};





/********************************************************************/
// A binning functor for accumulating a 2D histogram of angles and positions of molecules
// The functor can operate on a molecule so long as the molecule has a method for returning its molecular axis (vecr)
/*
	 template <class T>
	 class atomic_position_binner : public std::binary_function
	 public:

// The functor's operation for binning of atomic positions
void operator() (AtomPtr atom);

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
*/


/********************************************************************/
/********************************************************************/
/********************************************************************/

/*

	 template<class T>
	 class DensityAnalysis : public AnalysisSet<T> {

	 public:
	 typedef T system_t;

	 DensityAnalysis (std::string desc, std::string fn) : XYZAnalysisSet (desc,fn) { }
	 DensityAnalysis () :
	 AnalysisSet<system_t>() { return; }

// the utility functor for getting all the data accumulated
DensityBinner<Atom *, T> binner;

void Setup () { 
// Load all the atoms in the system 
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
*/

#endif
