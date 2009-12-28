#ifndef POS_ANG_HISTOGRAM_H_
#define POS_ANG_HISTOGRAM_H_

#include "../utility.h"
#include "../analysis.h"

//typedef vector<int> Histogram_1D;
//typedef vector<Histogram_1D> Histogram_2D;

/********************************************************************/
// A binning functor for accumulating a 2D histogram of angles and positions of molecules
// The functor can operate on a molecule so long as the molecule has a method for returning its molecular axis (vecr)
template <class T>
class APBinner {
  public:

	APBinner ();

	// The functor's operation for binning molecular orientation and position
	void operator() (T * mol) {
	  // calculate the bin for the cos(angle) formed between the molecular axis and the reference axis.
	  mol->SetOrderAxes();
	  //VecR molAxis = mol->MolecularAxis();
	  VecR molAxis = mol->Normal();
	  double angle = molAxis < Analyzer::ref_axis;
// flip equivalent angles (only good for when doing calculations on the molecular normal axis)
	  angle = (angle < 0.0) ? -angle : angle;
	  //int anglebin = Analyzer::AngleBin (molAxis < Analyzer::ref_axis);

	  // Calculate the bin for the molecule's position (based on center of mass)
	  mol->UpdateCenterOfMass();
	  double position = Analyzer::Position(mol->CenterOfMass());
	  //int positionbin = Analyzer::PositionBin (mol->CenterOfMass()[Analyzer::axis]);	

	  // perform the update/binning in the histogram
	  //_histogram.Shove (positionbin, anglebin);
	  //_histogram[positionbin][anglebin] += 1;
	  _histogram(position, angle);
	}

	// Return the histogram
	Histogram2D<double> Histogram () { return _histogram; }
	
	// Output data from a 2-dimensional histogram
	void APBinnerOutput (FILE * output) {

	  rewind (output);
	  std::pair<int,int> sizes = _histogram.size;

	  for (int i = 0; i < sizes.first; i++) {
		for (int j = 0; j < sizes.second; j++) {
		  fprintf (output, "%12d", _histogram.Element(i,j));
		}
		fprintf (output, "\n");
	  }
	  fflush (output);
	  return;
	}

	void PrintHistogram () {

	  std::pair<int,int> sizes = _histogram.size;

	  for (int i = 0; i < sizes.first; i++) {
		for (int j = 0; j < sizes.second; j++) {
			printf ("% 8d", _histogram.Element(i,j));
		}
	  }
	  fflush (stdout);
	}


  private:
	static Histogram2D<double> _histogram;		// The histogram for accumulating results
	static bool APBinner_instance;
};

template <class T> bool APBinner<T>::APBinner_instance = false;
template <class T> Histogram2D<double> APBinner<T>::_histogram (
	Histogram2D<double>::pair_t (150.0, 1.0),
	Histogram2D<double>::pair_t (0.1, 0.01),
	Histogram2D<double>::pair_t (-20.0, 0.0));
	//Histogram2D<double>::pair_t (WaterSystem<AmberSystem>::posmax, Analyzer::angmax),// maxima for position/angles
	//Histogram2D<double>::pair_t (Analyzer::posres, Analyzer::angres),	// resolutions
	//Histogram2D<double>::pair_t (WaterSystem<AmberSystem>::posmin, Analyzer::angmin));

// singleton-style constructor for the analysis functor
template <class T>
APBinner<T>::APBinner () {
  if (!APBinner<T>::APBinner_instance) {
	// make this the only time it gets set
	APBinner<T>::APBinner_instance = true;
  }
}


/********************************************************************/
class Angle_Position_Analyzer : public Analyzer {

  public:
	Angle_Position_Analyzer (WaterSystemParams& wsp) :
	  Analyzer(wsp) { }

	// the utility functor for getting all the data accumulated
	APBinner<Water> binner;

	void Setup () {
	  // For now let's look at waters only
	  FindWaters ();
	  return;
	}

	void Analysis () {
	  // bin each water's position and angle
	  FOR_EACH(int_wats, this->binner);
	  return;
	}

	// nothing at the moment
	void PostAnalysis () { return; }

	void DataOutput (const unsigned int timestep) {

	  if (!(timestep % (output_freq * 10)))
		this->binner.APBinnerOutput(output);
	  return;
	}

};
#endif
