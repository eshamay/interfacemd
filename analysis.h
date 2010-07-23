#ifndef ANALYSIS_H_
#define ANALYSIS_H_

#include "watersystem.h"
#include "utility.h"

template <class T>
class Analyzer : public WaterSystem<T> {

  protected:

    std::string output_filename;
    FILE * output;
    int	output_freq;

    void _OutputHeader () const;
    void _OutputStatus (const int timestep);
    void _CheckOutputFile ();

    void _EmptyFunction () const { return; } /* A simple empty function that does nothing to the system */


  public:
    Analyzer (const WaterSystemParams& params);
    virtual ~Analyzer ();


    virtual void SystemAnalysis ();

    static VecR ref_axis;

    // position boundaries and bin widths for gathering histogram data
    static double	posres;
    static int		posbins;
    static double 	angmin, angmax, angres;
    static int		angbins;
	int 			timestep;
    static int 		timesteps;
    static unsigned int restart;

    /*
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

    */

    // calculate a bin for a histogram
    static int Bin (const double value, const double min, const double res) {
      return (int)((value-min)/res);
    }

    static int PositionBin (const double position);
    static int PositionBin (const AtomPtr patom);
    static double Position (const AtomPtr patom);
    static double Position (const VecR& v);
    static double Position (const double d);

    static int AngleBin (const double angle) {
      return Bin (angle, angmin, angres);
    }

    void LoadNext ();


    Atom_ptr_vec& Atoms () { return WaterSystem<T>::int_atoms; } 
    Mol_ptr_vec& Molecules () { return WaterSystem<T>::int_mols; }
    Water_ptr_vec& Waters () { return WaterSystem<T>::int_wats; }

	// calculate the system's center of mass
	VecR CenterOfMass (const Mol_ptr_vec& mols) const;
	template <typename U> VecR CenterOfMass (const std::vector<U>& mols) const;

    // analysis loop functions
	virtual void Setup () = 0;
	virtual void Analysis () = 0;
	virtual void PostAnalysis () = 0;
	virtual void DataOutput () = 0;


	/*
	// check if a water molecule is above a certain location in the system
	class WaterAbovePosition : public std::binary_function <MolPtr, double, bool> {
	  public:
		bool operator () (const MolPtr m, const double cutoff) const {
		  return m->GetAtom("O")->Position()[WaterSystem<T>::axis] > cutoff;
		}
	};	// above position
	*/


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
Analyzer<T>::Analyzer (const WaterSystemParams& params)

: WaterSystem<T>(params),
  output_filename(params.output_filename), output(params.output),
  output_freq(params.output_freq),
  timestep (0)
{ 
  Analyzer<T>::posres = params.posres;
  Analyzer<T>::posbins = int((params.posmax - params.posmin)/params.posres);

  Analyzer<T>::angmin = params.angmin;
  Analyzer<T>::angmax = params.angmax;
  Analyzer<T>::angres = params.angres;
  Analyzer<T>::angbins = int((params.angmax - params.angmin)/params.angres);

  Analyzer<T>::timesteps = params.timesteps;
  Analyzer<T>::restart = params.restart;

  Analyzer<T>::ref_axis = params.ref_axis;

  this->_CheckOutputFile();
  this->_OutputHeader();
}

template <class T> 
void Analyzer<T>::_OutputHeader () const {

  printf ("Analysis Parameters:\n\tOutput Filename = \"%s\"\n\tScreen output frequency = 1/%d\n\n\tPosition extents for analysis:\n\t\tMin = % 8.3f\n\t\tMax = % 8.3f\n\t\tPosition Resolution = % 8.3f\n\n\tPrimary Axis = %d\nNumber of timesteps to be analyzed = %d\n",
	  output_filename.c_str(), output_freq, Analyzer<T>::posmin, Analyzer<T>::posmax, Analyzer<T>::posres, int(Analyzer<T>::axis), Analyzer<T>::timesteps);

#ifdef AVG
  printf ("\n\nThe analysis is averaging about the two interfaces located as:\n\tLow  = % 8.3f\n\tHigh = % 8.3f\n\n", int_low, int_high);
#endif
  return;
}

template <class T> 
Analyzer<T>::~Analyzer () {
  delete this->sys;
  return;
}

template <class T> 
void Analyzer<T>::_CheckOutputFile () {

  if (output == (FILE *)NULL) {
	printf ("WaterSystem::WaterSystem (argc, argv) - couldn't open the data output file!\n");
	exit(1);
  }

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
void Analyzer<XYZSystem>::LoadNext () {
  this->sys->LoadNext();
  this->LoadAll();
  return;
}

template <>
void Analyzer<AmberSystem>::LoadNext () {
  this->sys->LoadNext();
  return;
}



// A routine that performs some type of typical/generic analysis of a system
  template <class T> 
void Analyzer<T>::SystemAnalysis ()
{
  // do some initial setup
  this->Setup();

  // start the analysis - run through each timestep
  for (timestep = 0; timestep < timesteps; timestep++) {

	try {
	  // Perform the main loop analysis that works on every timestep of the simulation
	  this->Analysis ();
	} catch (std::exception& ex) {
	  std::cout << "Caught an exception during the system analysis at timestep " << timestep << "." << std::endl;
	  throw;
	}

	// output the status of the analysis (to the screen or somewhere useful)
	this->_OutputStatus (timestep);
	// Output the actual data being collected to a file or something for processing later
	if (!(timestep % (output_freq * 10)) && timestep)
	  this->DataOutput();


	try {
	  // load the next timestep
	  this->LoadNext();
	} catch (std::exception& ex) {
	  std::cout << "Caught an exception while loading the next timestep after step " << timestep << "." << std::endl;
	  throw;
	}
  }

  // do one final data output to push out the finalized data set
  DataOutput();

  // do a little work after the main analysis loop (normalization of a histogram? etc.)
  PostAnalysis ();
  return;
}	// system analysis

template <class T> 
int Analyzer<T>::PositionBin (const double position) {

  // check for flipping due to periodic boundaries
  double pos = position;
  if (pos < WaterSystem<T>::wsp.pbcflip) 
	pos += MDSystem::Dimensions()[WaterSystem<T>::wsp.axis];
  return (Bin (pos, WaterSystem<T>::wsp.posmin, posres));

}

// Find the position a particular atom fits into for a histogram
template <class T> 
int Analyzer<T>::PositionBin (const AtomPtr patom) {

  double position = Analyzer<T>::Position(patom);
  //VecR r = patom->Position();
  //double position = r[axis];
  //if (position < pbcflip) position += MDSystem::Dimensions()[axis];

#ifdef AVG
  // here the bin will be selected based on the distance to a given interface. Negative distances are inside the water phase, positive are in the CCl4
  double distance = (position > middle) ? position - int_high : int_low - position;
  int bin = Bin (distance, WaterSystem<T>::wsp.posmin, posres);
#else
  //if (position < START or position > END) continue;		// only bin stuff within the bounds that have been set up
  int bin = Bin (position, WaterSystem<T>::wsp.posmin, posres);
#endif

  return bin;
}

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
template <class U>
VecR Analyzer<T>::CenterOfMass (const std::vector<U>& mols) const 
{
  double mass = 0.0;
  VecR com;

  typedef typename std::vector<U>::const_iterator u_it;

  for (u_it it = mols.begin(); it != mols.end(); it++) {
	for (Atom_it jt = (*it)->begin(); jt != (*it)->end(); jt++) {
	  mass += (*jt)->Mass();
	  com += (*jt)->Position() * (*jt)->Mass();
	}
  }
  com /= mass;
  return com;
}



/*
// Create a histogram of the angles formed by an axis of a molecule's given reference axis.
// The supplied axis function should return the molecular axis vector of the molecule.
template <class T, typename molecule_t>
vector<int> Analyzer<T>::Molecular_Axis_Orientation_Histogram (
string molName,
VecR (molecule_t::*axisFn)() // This axisFn must return the molecular axis of interest of the molecule
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
 */

/*
   template <typename molecule_t>
   vector<int> Analyzer<T>::Molecular_Axis_Orientation_Position_Histogram (
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

/*
// creates a 2D histogram showing the shape of the plane formed by a given atom of a molecule averaged over a simulation.
// The histogram will show either the population density at a particular in-plane point giving the population as a function of the two coordinates...
// or it can show the average location (in the surface-normal direction) of atoms in the plane.
template <class T, typename molecule_t>
vector< vector<double> > Analyzer<T>::Interface_Location_Histogram (
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
 */

#endif
