#ifndef ANALYSIS_H_
#define ANALYSIS_H_

#include "watersystem.h"

template<class T>
class Analyzer : public WaterSystem<AmberSystem> {

  public:
	Analyzer (T& analysis_params, WaterSystemParams& params);

	void SystemAnalysis ();

	// create a histogram of the angle between a given molecular axis vector (determined by the axisFunc) and the system's ref_axis. The molecule is chosen by the residue name. The molecules must themselves have the functions for determining the molecular axis vector.
	template <typename molecule_t>
	vector<int> Molecular_Axis_Orientation_Histogram (string name) const;

	int AngleBin (const double angle) const;
	int PositionBin (const double position) const;
	int PositionBin (const Atom * patom) const;

  protected:
	void _OutputHeader () const;
	void _OutputStatus (const int timestep) const;

	// parameters for this particular analysis
	T _ap;

	// function to perform some initial setup before the main analysis loop
	virtual void Setup () = 0;
	// the main analysis function - this gets run on every timestep
	virtual void Analysis () = 0;
	// something to do after the loop (normalization, etc.) - done after the last timestep
	virtual void PostAnalysis () = 0;
	// A function that defines how data is output from the program
	virtual void DataOutput (const int timestep) = 0;
};

 
template <class T>
Analyzer<T>::Analyzer (T& ap, WaterSystemParams& wsp)
: WaterSystem<AmberSystem>(wsp), _ap(ap)
{ 
  this->sys = new AmberSystem("prmtop", "mdcrd", "mdvel");
  this->_OutputHeader();
}

template <class T>
void Analyzer<T>::_OutputHeader () const {

  printf ("Analysis Parameters:\n\tOutput Filename = \"%s\"\n\tScreen output frequency = 1/%d\n\n\tPosition extents for analysis:\n\t\tMin = % 8.3f\n\t\tMax = % 8.3f\n\t\tPosition Resolution = % 8.3f\n\n\tPrimary Axis = %d\nNumber of timesteps to be analyzed = %d\n",
	  output_filename.c_str(), output_freq, posmin, posmax, posres, int(axis), timesteps);

#ifdef AVG
  printf ("\n\nThe analysis is averaging about the two interfaces located as:\n\tLow  = % 8.3f\n\tHigh = % 8.3f\n\n", int_low, int_high);
#endif
  return;
}

template <class T>
void Analyzer<T>::_OutputStatus (const int timestep) const
{
  if (!(timestep % (this->output_freq * 10)))
	cout << endl << timestep << "/" << this->timesteps << " ) ";
  if (!(timestep % this->output_freq))
	cout << "*";

  fflush (stdout);
  return;
}

// A routine that performs some type of typical/generic analysis of a system
template <class T>
void Analyzer<T>::SystemAnalysis ()
{
  // do some initial setup
  Setup ();

  unsigned int timestep = 0;
  // start the analysis - run through each timestep
  for (timestep = 0; timestep < this->timesteps; timestep++) {

	// Perform the main loop analysis that works on every timestep of the simulation
	Analysis ();
	
	// load the next timestep
	this->sys->LoadNext();

	// output the status of the analysis (to the screen or somewhere useful)
	_OutputStatus (timestep);
	// Output the actual data being collected to a file or something for processing later
	DataOutput (timestep);
  }

  // do a little work after the main analysis loop (normalization of a histogram? etc.)
  PostAnalysis ();

  // do one final data output to push out the finalized data set
  DataOutput (timestep);

  return;
}

// calculate the histogram bin to be used given a particular angle
template <class T>
int Analyzer<T>::AngleBin (const double angle) const {
  int anglebin = int ((angle - angmin)/angres);
  return anglebin;
}

template <class T>
int Analyzer<T>::PositionBin (const double position) const {
  int posbin = int ((position - posmin)/posres);
  return posbin;
}

// Find the position a particular atom fits into for a histogram
template <class T>
int Analyzer<T>::PositionBin (const Atom * patom) const {

  VecR r = patom->Position();
  double position = r[axis];
  if (position < pbcflip) position += Atom::Size()[axis];

#ifdef AVG
  // here the bin will be selected based on the distance to a given interface. Negative distances are inside the water phase, positive are in the CCl4
  double distance = (position > middle) ? position - int_high : int_low - position;
  int bin = PositionBin (distance);
#else
  //if (position < START or position > END) continue;		// only bin stuff within the bounds that have been set up
  int bin = PositionBin (position);
#endif

  return bin;
}

// Create a histogram of the angles formed by an axis of a molecule's given reference axis. The supplied axis function should return the molecular axis vector of the molecule.
template <class T>
template <typename molecule_t>
vector<int> Analyzer<T>::Molecular_Axis_Orientation_Histogram (string molName)
const {

  // set up the histogram for output
  vector<int> histo (angbins, 0);

  molecule_t * pmol;
  // Run an analysis on all the carbon-chain molecules in the system to find their orientations over the course of a simulation with respect to a given axis.
  RUN (sys->Molecules()) {
	pmol = static_cast<molecule_t *>(sys->Molecules(i));
	// Search for the molNamed molecules in the system
	if (pmol->Name() != molName) continue;

	// find the particular molecular-axis vector
	VecR molAxis = pmol->Vector_CoM_To_End();
	double angle = (molAxis < ref_axis);	// calculates the cos(angle) between the two axes
	int anglebin = AngleBin (angle);
	histo[anglebin]++;
  }

  return (histo);
}
#endif
