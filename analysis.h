#ifndef ANALYSIS_H_
#define ANALYSIS_H_

#include "watersystem.h"

class Analyzer : public WaterSystem<AmberSystem> {

  protected:

	std::string output_filename;
	FILE * output;
	int	output_freq;

	void _OutputHeader () const;
	void _OutputStatus (const int timestep) const;
	void _CheckOutputFile ();

	void _EmptyFunction () const { return; } /* A simple empty function that does nothing to the system */



  public:
	Analyzer (WaterSystemParams& params);
	virtual ~Analyzer ();

	virtual void SystemAnalysis ();

	static VecR ref_axis;

	// position boundaries and bin widths for gathering histogram data
	static double	posres;
	static int		posbins;
	static double 	angmin, angmax, angres;
	static int		angbins;
	static unsigned int timesteps;
	static unsigned int restart;

	// create a histogram of the angle between a given molecular axis vector (determined by the axisFunc) and the system's ref_axis. The molecule is chosen by the residue name. The molecules must themselves have the functions for determining the molecular axis vector.
	template <typename molecule_t>
	  vector<int> Molecular_Axis_Orientation_Histogram (
		  string moleculeName,
		  VecR (molecule_t::*axisFunction)()
		  );

	// Creates a 2-D histogram of position in a slab vs. molecular orientation (calculated by the axis returned by the axisFunction). The vector returned can be accessed by vector[position][angle].
	// The true angles are not computed, but rather the cosines of the angles formed between the molecular axis and the system-reference axis.
	template <typename molecule_t>
	  vector< vector<int> > Molecular_Axis_Orientation_Position_Histogram (
		  string moleculeName,
		  VecR (molecule_t::*axisFunction)()
		  );

	template <typename molecule_t>
	  vector< vector<double> > Interface_Location_Histogram (
		  const coord d1, const coord d2,
		  const string mol_name,
		  const string atom_name,
		  vector< vector<double> >& histogram,
		  vector< vector<int> >& density);

	// calculate a bin for a histogram
	static int Bin (const double value, const double min, const double res) {
	  return (int)((value-min)/res);
	}

	static int PositionBin (const double position);
	static int PositionBin (const Atom * patom) ;

	static int AngleBin (const double angle) {
	  return Bin (angle, angmin, angres);
	}

	// analysis loop functions
	virtual void Setup () { return; }
	virtual void Analysis () { return; }
	virtual void PostAnalysis () { return; }
	virtual void DataOutput (const unsigned int timestep) { return; }

};


double	Analyzer::posres;
int		Analyzer::posbins;

double 	Analyzer::angmin;
double	Analyzer::angmax;
double	Analyzer::angres;
int		Analyzer::angbins;

unsigned int Analyzer::timesteps;
unsigned int Analyzer::restart;

VecR	Analyzer::ref_axis;


Analyzer::Analyzer (WaterSystemParams& params)

: WaterSystem<AmberSystem>(params),
  output_filename(params.output_filename), output(params.output),
  output_freq(params.output_freq)

{ 
  Analyzer::posres = params.posres;
  Analyzer::posbins = int((params.posmax - params.posmin)/params.posres);

  Analyzer::angmin = params.angmin;
  Analyzer::angmax = params.angmax;
  Analyzer::angres = params.angres;
  Analyzer::angbins = int((params.angmax - params.angmin)/params.angres);

  Analyzer::timesteps = params.timesteps;
  Analyzer::restart = params.restart;

  Analyzer::ref_axis = params.ref_axis;

  this->_CheckOutputFile();
  this->sys = new AmberSystem("prmtop", "mdcrd", "mdvel");
  this->_OutputHeader();
}

void Analyzer::_OutputHeader () const {

  printf ("Analysis Parameters:\n\tOutput Filename = \"%s\"\n\tScreen output frequency = 1/%d\n\n\tPosition extents for analysis:\n\t\tMin = % 8.3f\n\t\tMax = % 8.3f\n\t\tPosition Resolution = % 8.3f\n\n\tPrimary Axis = %d\nNumber of timesteps to be analyzed = %d\n",
	  output_filename.c_str(), output_freq, posmin, posmax, Analyzer::posres, int(Analyzer::axis), Analyzer::timesteps);

#ifdef AVG
  printf ("\n\nThe analysis is averaging about the two interfaces located as:\n\tLow  = % 8.3f\n\tHigh = % 8.3f\n\n", int_low, int_high);
#endif
  return;
}

Analyzer::~Analyzer () {
  fclose (output);
  return;
}

void Analyzer::_CheckOutputFile () {

  if (output == (FILE *)NULL) {
	printf ("WaterSystem::WaterSystem (argc, argv) - couldn't open the data output file!\n");
	exit(1);
  }

  return;
}

void Analyzer::_OutputStatus (const int timestep) const
{
  if (!(timestep % (this->output_freq * 10)))
	cout << endl << timestep << "/" << this->timesteps << " ) ";
  if (!(timestep % this->output_freq))
	cout << "*";

  fflush (stdout);
  return;
}

// A routine that performs some type of typical/generic analysis of a system
void Analyzer::SystemAnalysis ()
{
  // do some initial setup
  this->Setup();

  unsigned int timestep = 0;
  // start the analysis - run through each timestep
  for (timestep = 0; timestep < this->timesteps; timestep++) {

	// Perform the main loop analysis that works on every timestep of the simulation
	this->Analysis ();

	// load the next timestep
	this->sys->LoadNext();

	// output the status of the analysis (to the screen or somewhere useful)
	this->_OutputStatus (timestep);
	// Output the actual data being collected to a file or something for processing later
	this->DataOutput(timestep);
  }

  // do a little work after the main analysis loop (normalization of a histogram? etc.)
  this->PostAnalysis ();

  // do one final data output to push out the finalized data set
  this->DataOutput(timestep);

  return;
}

int Analyzer::PositionBin (const double position) {

  // check for flipping due to periodic boundaries
  double pos = position;
  if (pos < pbcflip) 
	pos += MDSystem::Dimensions()[axis];
  return (Bin (pos, posmin, posres));

}

// Find the position a particular atom fits into for a histogram
int Analyzer::PositionBin (const Atom * patom) {

  VecR r = patom->Position();
  double position = r[axis];
  if (position < pbcflip) position += MDSystem::Dimensions()[axis];

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
template <typename molecule_t>
vector<int> Analyzer::Molecular_Axis_Orientation_Histogram (
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
	const double angle = (molAxis < ref_axis);	// calculates the cos(angle) between the two axes
	int anglebin = this->Bin (angle, angmin, angres);
	histo[anglebin]++;
  }

  return (histo);
}

/*
template <typename molecule_t>
vector<int> Analyzer::Molecular_Axis_Orientation_Position_Histogram (
	string molName,
	VecR (molecule_t::*axisFunction)() 	// This axisFn must return the molecular axis of interest of the molecule
	)
{

  // set up a 2-dimensional histogram for output - histo[position-bin][angle-bin]
  vector< vector<int> > histo (posbins, vector<int> (angbins, 0));

  molecule_t * pmol;
  // Each molecule of interest is asked for its:
  RUN (int_mols) {
	pmol = static_cast<molecule_t *>(int_mols[i]);

	// position
	int posbin = PositionBin (pmol->GetAtom(positioningAtom));

	// and angle between the molecular axis of interest and the reference axis
	VecR molAxis = (pmol->*axisFn)();
	const double angle = (molAxis < ref_axis);	// calculates the cos(angle) between the two axes
	int anglebin = this->Bin (angle, angmin, angres);

	histo[posbin][anglebin]++;	// Binning into the histogram, ftw
  }

  return (histo);
}
*/

// creates a 2D histogram showing the shape of the plane formed by a given atom of a molecule averaged over a simulation.
// The histogram will show either the population density at a particular in-plane point giving the population as a function of the two coordinates...
// or it can show the average location (in the surface-normal direction) of atoms in the plane.
template <typename molecule_t>
vector< vector<double> > Analyzer::Interface_Location_Histogram (
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
  double d1_max = MDSystem::Dimensions()[d1] + 10.0;	
  // the size of the system in the d2 direction
  double d2_max = MDSystem::Dimensions()[d2] + 10.0;
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
