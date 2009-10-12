#include "analysis.h"

Analyzer::Analyzer (
	void * analysis_params,
	WaterSystemParams& params
	)
: WaterSystem<AmberSystem>(params), _ap(analysis_params)
{ 
	this->sys = new AmberSystem("prmtop", "mdcrd", "mdvel");
	this->OutputHeader();
}

void Analyzer::_OutputHeader () const {

	printf ("Analysis Parameters:\n\tOutput Filename = \"%s\"\n\tScreen output frequency = 1/%d\n\n\tPosition extents for analysis:\n\t\tMin = % 8.3f\n\t\tMax = % 8.3f\n\t\tPosition Resolution = % 8.3f\n\n\tPrimary Axis = %d\nNumber of timesteps to be analyzed = %d\n",
			output.c_str(), output_freq, posmin, posmax, posres, int(axis), timesteps);

	if (params.avg) {
		printf ("\n\nThe analysis is averaging about the two interfaces located as:\n\tLow  = % 8.3f\n\tHigh = % 8.3f\n\n", int_low, int_high);
	}

	return;
}

void Analyzer::_OutputStatus (const int timestep)
{
  if (!(timestep % (this->output_freq * 10)))
	cout << endl << timestep << "/" << this->timesteps << " ) ";
  if (!(timestep % this->output_freq))
	cout << "*";

  fflush (stdout);
  return;
}

// A routine that performs some type of typical/generic analysis of a system
void SystemAnalysis ()
{
  // do some initial setup
  Setup ();

  // start the analysis - run through each timestep
  for (timestep = 0; timestep < this->timesteps; timestep++) {

	// Perform the analysis
	Analysis ();

	this->sys->LoadNext();

	// output the status of the analysis (to the screen or somewhere useful)
	_OutputStatus (timestep);
	// Output the actual data being collected to a file or something
	DataOutput (timestep);
  }

  // do a little work after the main analysis loop
  PostAnalysis ();

  // do one final data output to push out the finalized data set
  DataOutput (timestep);

  return;
}

// calculate the histogram bin to be used given a particular angle
int Analyzer::AngleBin (const double angle) {
  int anglebin = int ((angle - anglemin)/angleres);
  return anglebin;
}

int Analyzer::PositionBin (const double position) {
  int posbin = int ((pos - posmin)/posres);
  return posbin;
}

// Find the position a particular atom fits into for a histogram
int Analyzer::PositionBin (const Atom * patom) {

  VecR r = patom->Position();
  double position = r[ref_axis];
  if (position < pbcflip) position += Atom::Size()[ref_axis];

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
vector<int> Analyzer::Molecular_Axis_Orientation_Histogram (
	const string name,
	VecR (*axisFunc)()
	)
{

  // set up the histogram for output
  vector<int> histo (angbins, 0);

  Molecule * pmol;
  // Run an analysis on all the carbon-chain molecules in the system to find their orientations over the course of a simulation with respect to a given axis.
  RUN (sys->Molecules()) {
	pmol = sys->Molecules(i);
	// Search for the named molecules in the system
	if (pmol->Name() != name) continue;

	// find the particular molecular-axis vector
	VecR molAxis = pmol->axisFunc();
	double angle = (molAxis < ref_axis);

	int anglebin = AngleBin (angle);
	histo[anglebin]++;
  }

  return (histo);
}
