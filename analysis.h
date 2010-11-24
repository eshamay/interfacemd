#ifndef ANALYSIS_H_
#define ANALYSIS_H_

#include "watersystem.h"
#include "utility.h"
#include "patterns.h"
#include "dataoutput.h"



// An analysis that will be performed on a system by an analyzer
template <typename T>
class AnalysisSet {
	public:
		typedef T system_t;

		virtual ~AnalysisSet () { }
		AnalysisSet (std::string desc, std::string fn) : description (desc), filename(fn) { }

		// default setup
		virtual void Setup (system_t& t) {
			rewind(t.Output());
			t.LoadAll();
			return;
		}

		// each analyzer has to have an analysis function to do some number crunching
		virtual void Analysis (system_t&) = 0;
		// normally this can be done in the analysis section, but just for style we can have something different defined here
		virtual void DataOutput (system_t&) { }
		virtual void PostAnalysis (system_t&) { }

		std::string& Description () { return description; }
		std::string& Filename () { return filename; }

	protected:

		std::string description;	// describes the analysis that is performed
		std::string filename;		// filename to use for data output

};  // class AnalysisSet



template <class T>
class Analyzer : public WaterSystem<T>, public patterns::observer::observable {

	protected:

		std::string output_filename;
		FILE * output;
		int	output_freq;

		void _OutputHeader () const;
		void _OutputStatus ();
		md_analysis::StarStatusBarUpdater	status_updater;

	public:
		Analyzer (const std::string = std::string("system.cfg"));
		virtual ~Analyzer ();

		typedef AnalysisSet<Analyzer<T> > analysis_t;
		void SystemAnalysis (analysis_t&);

		// position boundaries and bin widths for gathering histogram data
		static double	posres;
		static int		posbins;
		static double angmin, angmax, angres;
		static int		angbins;
		int						timestep;
		int						Timestep () const { return timestep; }
		static int 		timesteps;
		static unsigned int restart;

		static double Position (const AtomPtr);
		static double Position (const VecR&);
		static double Position (const double);

		void LoadNext ();

		FILE * Output () { return output; }
		void OpenDataOutputFile (analysis_t&);

		Atom_ptr_vec& Atoms () { return WaterSystem<T>::int_atoms; } 
		Mol_ptr_vec& Molecules () { return WaterSystem<T>::int_mols; }
		Water_ptr_vec& Waters () { return WaterSystem<T>::int_wats; }

		// calculate the system's center of mass
		template <typename Iter> static VecR CenterOfMass (Iter first, Iter last);

		//! Predicate for sorting a container of molecules based on position along the main axis of the system, and using a specific element type to determine molecular position. i.e. sort a container of waters based on the O position, or sort a container of NO3s based on the N position, etc.
		class molecule_position_pred; 		
		class molecule_reference_distance_pred; 		
		class atomic_reference_distance_pred;

};	// Analyzer


template <class T> double	Analyzer<T>::posres;
template <class T> int		Analyzer<T>::posbins;

template <class T> double 	Analyzer<T>::angmin;
template <class T> double	Analyzer<T>::angmax;
template <class T> double	Analyzer<T>::angres;
template <class T> int		Analyzer<T>::angbins;

template <class T> int 		Analyzer<T>::timesteps;
template <class T> unsigned int Analyzer<T>::restart;

template <class T> 
Analyzer<T>::Analyzer (const std::string ConfigurationFilename) : 
	WaterSystem<T>(ConfigurationFilename),
	output_filename(""), output((FILE *)NULL),
	output_freq(WaterSystem<T>::SystemParameterLookup("analysis.output-frequency")),
	timestep (0)
{ 
	Analyzer<T>::posres = WaterSystem<T>::SystemParameterLookup("analysis.resolution.position");
	Analyzer<T>::posbins = int((WaterSystem<T>::posmax - WaterSystem<T>::posmin)/posres);

	Analyzer<T>::angmin = WaterSystem<T>::SystemParameterLookup("analysis.angle-range")[0];
	Analyzer<T>::angmax = WaterSystem<T>::SystemParameterLookup("analysis.angle-range")[1];
	Analyzer<T>::angres = WaterSystem<T>::SystemParameterLookup("analysis.resolution.angle");
	Analyzer<T>::angbins = int((angmax - angmin)/angres);

	Analyzer<T>::timesteps = WaterSystem<T>::SystemParameterLookup("system.timesteps");
	Analyzer<T>::restart = WaterSystem<T>::SystemParameterLookup("analysis.restart-time");

	status_updater.Set (output_freq, timesteps, 0);
	this->registerObserver(&status_updater);

	this->_OutputHeader();
} // Analyzer ctor



template <class T> 
void Analyzer<T>::_OutputHeader () const {

	printf ("Analysis Parameters:\n\tScreen output frequency = 1/%d\n\n\tPosition extents for analysis:\n\t\tMin = % 8.3f\n\t\tMax = % 8.3f\n\t\tPosition Resolution = % 8.3f\n\n\tPrimary Axis = %d\nNumber of timesteps to be analyzed = %d\n",
			output_freq, Analyzer<T>::posmin, Analyzer<T>::posmax, Analyzer<T>::posres, int(Analyzer<T>::axis), Analyzer<T>::timesteps);

#ifdef AVG
	printf ("\n\nThe analysis is averaging about the two interfaces located as:\n\tLow  = % 8.3f\n\tHigh = % 8.3f\n\n", int_low, int_high);
#endif
	return;
}

template <class T> 
Analyzer<T>::~Analyzer () {

	delete this->sys;
	if (output != (FILE *)NULL)
		fclose(output);

	return;
}

template <class T> 
void Analyzer<T>::OpenDataOutputFile (analysis_t& an) {

	output = (FILE *)NULL;
	output = fopen(an.Filename().c_str(), "w");

	if (output == (FILE *)NULL) {
		printf ("Analyzer<T>::_CheckOutputFile() - couldn't open the data output file, \"%s\", given in the analysis set!\n", an.Filename().c_str());
		exit(1);
	}

	printf ("\nOutputting data to \"%s\"\n", an.Filename().c_str());

	return;
}

	template <class T> 
void Analyzer<T>::_OutputStatus ()
{
	this->notifyObservers ();
	/*
	if (!(timestep % (this->output_freq * 10)))
		std::cout << std::endl << timestep << "/" << this->timesteps << " ) ";
	if (!(timestep % this->output_freq)) {
		std::cout << "*";
	}

	fflush (stdout);
	*/
	return;
}


template <>
extern void Analyzer<XYZSystem>::LoadNext () {
	this->sys->LoadNext();
	this->LoadAll();
	return;
}

template <>
extern void Analyzer<AmberSystem>::LoadNext () {
	this->sys->LoadNext();
	return;
}

template <>
extern void Analyzer<gromacs::GMXSystem<gromacs::TRRFile> >::LoadNext () {
	this->sys->LoadNext();
	return;
}

template <>
extern void Analyzer<gromacs::GMXSystem<gromacs::XTCFile> >::LoadNext () {
	this->sys->LoadNext();
	return;
}

template <class T>
void Analyzer<T>::SystemAnalysis (analysis_t& an) {
	// Open a file for data output
	this->OpenDataOutputFile (an);
	// do some initial setup
	an.Setup(*this);

	// start the analysis - run through each timestep
	for (timestep = 0; timestep < timesteps; timestep++) {

		try {
			// Perform the main loop analysis that works on every timestep of the simulation
			an.Analysis (*this);
		} catch (std::exception& ex) {
			std::cout << "Caught an exception during the system analysis at timestep " << timestep << "." << std::endl;
			throw;
		}

		// output the status of the analysis (to the screen or somewhere useful)
		this->_OutputStatus ();
		// Output the actual data being collected to a file or something for processing later
		if (!(timestep % (output_freq * 10)) && timestep)
			an.DataOutput(*this);


		try {
			// load the next timestep
			this->LoadNext();
		} catch (std::exception& ex) {
			throw;
		}
	}

	// do one final data output to push out the finalized data set
	an.DataOutput(*this);

	// do a little work after the main analysis loop (normalization of a histogram? etc.)
	an.PostAnalysis (*this);
	return;
} // System Analysis w/ analysis set



/* Find the periodic-boundary-satistfying location of an atom, vector, or raw coordinate along the reference axis */
template <class T> 
double Analyzer<T>::Position (const AtomPtr patom) {
	//return Analyzer<T>::Position(patom->Position());
	return WaterSystem<T>::AxisPosition(patom);
}

template <class T> 
double Analyzer<T>::Position (const VecR& v) {
	double position = v[WaterSystem<T>::axis];
	return Analyzer<T>::Position(position);
}

template <class T> 
double Analyzer<T>::Position (const double d) {
	double pos = d;
	if (pos < WaterSystem<T>::pbcflip) pos += MDSystem::Dimensions()[WaterSystem<T>::axis];
	return pos;
}


template <class T>
	template <class Iter>	// Has to be iterators to a container of molecules
VecR Analyzer<T>::CenterOfMass (Iter first, Iter last)
{
	double mass = 0.0;
	VecR com;
	com.setZero();

	typedef typename std::iterator_traits<Iter>::value_type val_t;

	for (Iter it = first; it != last; it++) {
		for (Atom_it jt = (*it)->begin(); jt != (*it)->end(); jt++) {
			mass += (*jt)->Mass();
			com += (*jt)->Position() * (*jt)->Mass();
		}
	}
	com /= mass;
	return com;
}

//! Predicate for sorting molecules based on their positions along the system reference axis. The position of the element supplied (elmt) is used. e.g. if elmt = Atom::O, then the first oxygen of the molecule will be used
template<class T>
class Analyzer<T>::molecule_position_pred : public std::binary_function <Molecule*,Molecule*,bool> {
	private:
		Atom::Element_t _elmt;	// determines the element in a molecule to use for position comparison
	public:
		//! upon instantiation, the element to be used for specifying molecular position is provided
		molecule_position_pred (const Atom::Element_t elmt) : _elmt(elmt) { }

		bool operator()(const Molecule* left, const Molecule* right) const {
			AtomPtr left_o = left->GetAtom(_elmt);
			AtomPtr right_o = right->GetAtom(_elmt);
			double left_pos = Analyzer<T>::Position(left_o);
			double right_pos = Analyzer<T>::Position(right_o);

			return left_pos < right_pos;
		}
};


//! Calculates the distance between two molecules using their given reference points (as defined in the Molecule::ReferencePoint
template<class T>
class Analyzer<T>::molecule_reference_distance_pred : public std::binary_function <MolPtr,MolPtr,bool> {
	private:
		MolPtr _reference_mol;	// the molecule that will act as the reference point for the comparison
	public:
		molecule_reference_distance_pred (const MolPtr refmol) : _reference_mol(refmol) { }
		// return the distance between the two molecules and the reference mol
		bool operator()(const MolPtr left, const MolPtr right) const {
			double left_dist = (left->ReferencePoint() - _reference_mol->ReferencePoint()).norm();
			double right_dist = (right->ReferencePoint() - _reference_mol->ReferencePoint()).norm();
			return left_dist < right_dist;
		}
};

// this predicate is used for distance calculations/sorting between atoms given a reference atom
template<class T>
class Analyzer<T>::atomic_reference_distance_pred : public std::binary_function <AtomPtr,AtomPtr,bool> {
	private:
		AtomPtr _refatom;	// the molecule that will act as the reference point for the comparison
	public:
		atomic_reference_distance_pred (const AtomPtr refatom) : _refatom (refatom) { }
		// return the distance between the two molecules and the reference mol
		bool operator()(const AtomPtr left, const AtomPtr right) const {
			double left_dist = (left->Position() - _refatom->Position()).norm();
			double right_dist = (right->Position() - _refatom->Position()).norm();
			return left_dist < right_dist;
		}
};

/***************** Analysis Sets specific to given MD systems ***************/
class XYZAnalysisSet : public AnalysisSet< Analyzer<XYZSystem> > { 
	public:
		typedef Analyzer<XYZSystem> system_t;
		XYZAnalysisSet (std::string desc, std::string fn) :
			AnalysisSet<system_t> (desc, fn) { }
		virtual ~XYZAnalysisSet () { }
};

class AmberAnalysisSet : public AnalysisSet< Analyzer<AmberSystem> > { 
	public:
		typedef Analyzer<AmberSystem> system_t;
		AmberAnalysisSet (std::string desc, std::string fn) :
			AnalysisSet< system_t > (desc, fn) { }
		virtual ~AmberAnalysisSet () { }
};

#endif
