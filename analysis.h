#ifndef ANALYSIS_H_
#define ANALYSIS_H_

#include "watersystem.h"

template<class T, typename U>
class Analyzer : public WaterSystem<AmberSystem> {

  public:
	Analyzer (U& analysis_params, WaterSystemParams& params);
	~Analyzer ();

	void SystemAnalysis (T * system);

	// create a histogram of the angle between a given molecular axis vector (determined by the axisFunc) and the system's ref_axis. The molecule is chosen by the residue name. The molecules must themselves have the functions for determining the molecular axis vector.
	template <typename molecule_t>
		vector<int> Molecular_Axis_Orientation_Histogram (
			string name,
			VecR (molecule_t::*axisFn)()
			);

	template <typename molecule_t>
	vector< vector<double> > Interface_Location_Histogram (
			const coord d1, const coord d2,
			const string mol_name,
			const string atom_name,
			vector< vector<double> >& histogram,
			vector< vector<int> >& density);

	int AngleBin (const double angle) const;
	int PositionBin (const double position) const;
	int PositionBin (const Atom * patom) const;
	// calculate a bin for a histogram
	int Bin (const double value, const double min, const double res) const;

	void Set_Setup 			(void (T::*FnPtr)()) { _setup = FnPtr; }
	void Set_Analysis 		(void (T::*FnPtr)()) { _analysis = FnPtr; }
	void Set_PostAnalysis 	(void (T::*FnPtr)()) { _post_analysis = FnPtr; }
	void Set_DataOutput 	(void (T::*FnPtr)(const int timestep)) 
	{ 
		_data_output = FnPtr; 
	}

  protected:

	// parameters for this particular analysis
	U _ap;

	std::string output_filename;
	FILE * output;

	VecR ref_axis;

	int	output_freq;

	// position boundaries and bin width
	double	posres;
	int		posbins;
	double 	angmin, angmax, angres;
	int		angbins;
	unsigned int timesteps;
	unsigned int restart;

	void _OutputHeader () const;
	void _OutputStatus (const int timestep) const;

	void _EmptyFunction () const { return; } /* A simple empty function that does nothing to the system */

	void CheckOutputFile ();

	// function pointers to perform the actual system analysis operations
	typedef void (T::*FnPtr)();
	typedef void (T::*DataOutputPtr)(const int timestep);
	FnPtr _setup, _analysis, _post_analysis;
	DataOutputPtr	_data_output;

};


template <class T, typename U>
Analyzer<T,U>::Analyzer (U& ap, WaterSystemParams& wsp)

: WaterSystem<AmberSystem>(wsp), _ap(ap),
	output_filename(wsp.output_filename), output(wsp.output),
	ref_axis(wsp.ref_axis),
	output_freq(wsp.output_freq),
	posres(wsp.posres), posbins(int((posmax - posmin)/posres)),
	angmin(wsp.angmin), angmax(wsp.angmax), angres(wsp.angres),
	angbins(int ((angmax - angmin)/angres) + 1), 
	timesteps(wsp.timesteps), restart(wsp.restart)

{ 
	this->CheckOutputFile();
	this->sys = new AmberSystem("prmtop", "mdcrd", "mdvel");
	this->_OutputHeader();
}

template <class T, typename U>
void Analyzer<T,U>::_OutputHeader () const {

  printf ("Analysis Parameters:\n\tOutput Filename = \"%s\"\n\tScreen output frequency = 1/%d\n\n\tPosition extents for analysis:\n\t\tMin = % 8.3f\n\t\tMax = % 8.3f\n\t\tPosition Resolution = % 8.3f\n\n\tPrimary Axis = %d\nNumber of timesteps to be analyzed = %d\n",
	  output_filename.c_str(), output_freq, posmin, posmax, posres, int(axis), timesteps);

#ifdef AVG
  printf ("\n\nThe analysis is averaging about the two interfaces located as:\n\tLow  = % 8.3f\n\tHigh = % 8.3f\n\n", int_low, int_high);
#endif
  return;
}

template <class T, typename U>
Analyzer<T,U>::~Analyzer () {
	fclose (output);
	return;
}

template <class T, typename U>
void Analyzer<T,U>::CheckOutputFile () {

	if (output == (FILE *)NULL) {
		printf ("WaterSystem::WaterSystem (argc, argv) - couldn't open the data output file!\n");
		exit(1);
	}

	return;
}

template <class T, typename U>
void Analyzer<T,U>::_OutputStatus (const int timestep) const
{
  if (!(timestep % (this->output_freq * 10)))
	cout << endl << timestep << "/" << this->timesteps << " ) ";
  if (!(timestep % this->output_freq))
	cout << "*";

  fflush (stdout);
  return;
}

// A routine that performs some type of typical/generic analysis of a system
template <class T, typename U>
void Analyzer<T,U>::SystemAnalysis (T * system)
{
  // do some initial setup
  (system->*_setup) ();

  unsigned int timestep = 0;
  // start the analysis - run through each timestep
  for (timestep = 0; timestep < this->timesteps; timestep++) {

	// Perform the main loop analysis that works on every timestep of the simulation
	(system->*_analysis) ();
	
	// load the next timestep
	this->sys->LoadNext();

	// output the status of the analysis (to the screen or somewhere useful)
	_OutputStatus (timestep);
	// Output the actual data being collected to a file or something for processing later
	(system->*_data_output) (timestep);
  }

  // do a little work after the main analysis loop (normalization of a histogram? etc.)
  (system->*_post_analysis) ();

  // do one final data output to push out the finalized data set
  (system->*_data_output) (timestep);

  return;
}

template <class T, typename U>
// given a value that needs to fit into a histogram, the bin is calculated based on the minimum bin-value of the histogram, and the bin-width, or resoluation.
int Analyzer<T,U>::Bin (const double value, const double min, const double res) 
	const
{
	int bin = (int)((value-min)/res);
	return bin;
}

// Find the position a particular atom fits into for a histogram
template <class T, typename U>
int Analyzer<T,U>::PositionBin (const Atom * patom) const {

  VecR r = patom->Position();
  double position = r[axis];
  if (position < pbcflip) position += Atom::Size()[axis];

#ifdef AVG
  // here the bin will be selected based on the distance to a given interface. Negative distances are inside the water phase, positive are in the CCl4
  double distance = (position > middle) ? position - int_high : int_low - position;
  int bin = Bin (distance, posmin, posres);
#else
  //if (position < START or position > END) continue;		// only bin stuff within the bounds that have been set up
  int bin = Bin (position, posmin, posres);
#endif

  return bin;
}

// Create a histogram of the angles formed by an axis of a molecule's given reference axis.
// The supplied axis function should return the molecular axis vector of the molecule.
template <class T, typename U>
template <typename molecule_t>
vector<int> Analyzer<T,U>::Molecular_Axis_Orientation_Histogram (
		string molName,
		VecR (molecule_t::*axisFn)() /* This axisFn must return the molecular axis of interest of the molecule */
		)
{

	// set up the histogram for output
	vector<int> histo (angbins, 0);

	molecule_t * pmol;
	// Run an analysis on all the carbon-chain molecules in the system to find their orientations over the course of a simulation with respect to a given axis.
	RUN (int_mols) {
		pmol = static_cast<molecule_t *>(int_mols[i]);

		// find the particular molecular-axis vector
		VecR molAxis = (pmol->*axisFn)();
		double angle = (molAxis < ref_axis);	// calculates the cos(angle) between the two axes
		int anglebin = AngleBin (angle);
		histo[anglebin]++;
	}

  return (histo);
}

// creates a 2D histogram showing the shape of the plane formed by a given atom of a molecule averaged over a simulation.
// The histogram will show either the population density at a particular in-plane point giving the population as a function of the two coordinates...
// or it can show the average location (in the surface-normal direction) of atoms in the plane.
template <class T, typename U>
template <typename molecule_t>
vector< vector<double> > Analyzer<T,U>::Interface_Location_Histogram (
								      // 2 directions (axes) parallel to the plane
								      const coord d1, const coord d2,
								      // name of the molecule to analyze
								      const string mol_name,
								      // Atom name to analyze
								      const string atom_name,
								      vector< vector<double> >& histogram,
								      vector< vector<int> >& density)
{
  double d_min = -5.0;
  // the size of the system in the d1 direction
  double d1_max = Atom::Size()[d1] + 10.0;	
  // the size of the system in the d2 direction
  double d2_max = Atom::Size()[d2] + 10.0;
  double d_res = 1.0;

  int numBins1 = (int)((d1_max - d_min)/d_res);
  int numBins2 = (int)((d2_max - d_min)/d_res);

  if (histogram.size() == 0)
    histogram.resize (numBins1, vector<double> (numBins2, 0));
  if (density.size() == 0)
    density.resize (numBins1, vector<int> (numBins2, 0));
	
  molecule_t * mol;
  Atom * atom;
  int d1_bin, d2_bin;
  RUN (int_mols) {
    mol = static_cast<molecule_t *>(int_mols[i]);
    atom = mol->GetAtom(atom_name);
		
    d1_bin = Bin (atom->Position()[d1], d_min, d_res);
    d2_bin = Bin (atom->Position()[d2], d_min, d_res);
    // bin the location of an atom at each point on the plane
    double ref_position = atom->Position()[axis];

    histogram[d1_bin][d2_bin] += ref_position;
    density[d1_bin][d2_bin]++;
  }

  return histogram;
}

#endif
