#ifndef ANALYSIS_H_
#define ANALYSIS_H_

#include "watersystem.h"
#include "utility.h"




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
class Analyzer : public WaterSystem<T> {

	protected:

		std::string output_filename;
		FILE * output;
		int	output_freq;

		void _OutputHeader () const;
		void _OutputStatus (const int);


	public:
		Analyzer (const std::string = std::string("system.cfg"));
		virtual ~Analyzer ();

		typedef AnalysisSet<Analyzer<T> > analysis_t;
		void SystemAnalysis (analysis_t&);

		static VecR ref_axis;

		// position boundaries and bin widths for gathering histogram data
		static double	posres;
		static int		posbins;
		static double 	angmin, angmax, angres;
		static int		angbins;
		int						timestep;
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

};	// Analyzer


template <class T> double	Analyzer<T>::posres;
template <class T> int		Analyzer<T>::posbins;

template <class T> double 	Analyzer<T>::angmin;
template <class T> double	Analyzer<T>::angmax;
template <class T> double	Analyzer<T>::angres;
template <class T> int		Analyzer<T>::angbins;

template <class T> int 		Analyzer<T>::timesteps;
template <class T> unsigned int Analyzer<T>::restart;

template <class T> VecR		Analyzer<T>::ref_axis;

template <class T> 
Analyzer<T>::Analyzer (const std::string ConfigurationFilename) : 
	WaterSystem<T>(ConfigurationFilename),
	output_filename(""), output((FILE *)NULL),
	output_freq(WaterSystem<T>::SystemParameters()->output_freq),
	timestep (0)
{ 
	Analyzer<T>::posres = WaterSystem<T>::SystemParameters()->posres;
	Analyzer<T>::posbins = int((WaterSystem<T>::SystemParameters()->posmax - WaterSystem<T>::SystemParameters()->posmin)/WaterSystem<T>::SystemParameters()->posres);

	Analyzer<T>::angmin = WaterSystem<T>::SystemParameters()->angmin;
	Analyzer<T>::angmax = WaterSystem<T>::SystemParameters()->angmax;
	Analyzer<T>::angres = WaterSystem<T>::SystemParameters()->angres;
	Analyzer<T>::angbins = int((WaterSystem<T>::SystemParameters()->angmax - WaterSystem<T>::SystemParameters()->angmin)/WaterSystem<T>::SystemParameters()->angres);

	Analyzer<T>::timesteps = WaterSystem<T>::SystemParameters()->timesteps;
	Analyzer<T>::restart = WaterSystem<T>::SystemParameters()->restart;

	Analyzer<T>::ref_axis = WaterSystem<T>::SystemParameters()->ref_axis;

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
void Analyzer<T>::_OutputStatus (const int timestep)
{
	if (!(timestep % (this->output_freq * 10)))
		std::cout << std::endl << timestep << "/" << this->timesteps << " ) ";
	if (!(timestep % this->output_freq)) {
		std::cout << "*";
	}

	fflush (stdout);
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
		this->_OutputStatus (timestep);
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
	return Analyzer<T>::Position(patom->Position());
}

template <class T> 
double Analyzer<T>::Position (const VecR& v) {
	double position = v[WaterSystem<T>::wsp.axis];
	return Analyzer<T>::Position(position);
}

template <class T> 
double Analyzer<T>::Position (const double d) {
	double pos = d;
	if (pos < WaterSystem<T>::wsp.pbcflip) pos += MDSystem::Dimensions()[WaterSystem<T>::wsp.axis];
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
