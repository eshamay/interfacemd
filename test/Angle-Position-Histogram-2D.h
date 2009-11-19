#ifndef POS_ANG_HISTOGRAM_H_
#define POS_ANG_HISTOGRAM_H_

#include "../utility.h"
#include "../analysis.h"

typedef vector<int> Histogram_1D;
typedef vector<Histogram_1D> Histogram_2D;

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
	  VecR molAxis = mol->Normal();
	  //VecR molAxis = mol->MolecularAxis();
	  int anglebin = Analyzer::AngleBin (molAxis < Analyzer::ref_axis);

	  // Calculate the bin for the molecule's position (based on center of mass)
	  mol->UpdateCenterOfMass();
	  int positionbin = Analyzer::PositionBin (mol->CenterOfMass()[Analyzer::axis]);	
	  // perform the update/binning in the histogram
	  _histogram[positionbin][anglebin] += 1;
	}

	// Return the histogram
	Histogram_2D Histogram () { return _histogram; }
	
	// Output data from a 2-dimensional histogram
	void HistogramOutput (FILE * output) {

	  rewind (output);

	  RUN(_histogram) {
		RUN2 (_histogram[i]) {
		  fprintf (output, "% 8d", _histogram[i][j]);
		}
		fprintf (output, "\n");
	  }
	  fflush (output);
	  return;
	}

	void PrintHistogram () {

	  RUN(_histogram) {
		RUN2 (_histogram[i]) {
		  if (_histogram[i][j] != 0)
			printf ("% 8d", _histogram[i][j]);
		}
	  }
	  fflush (stdout);
	}


  private:
	static Histogram_2D _histogram;		// The histogram for accumulating results
	static bool APBinner_instance;
};

template <class T> bool APBinner<T>::APBinner_instance = false;
template <class T> Histogram_2D APBinner<T>::_histogram;

// singleton-style constructor for the analysis functor
template <class T>
APBinner<T>::APBinner () {
  if (!APBinner<T>::APBinner_instance) {
	// create the histogram for binning the right amount of data
	APBinner<T>::_histogram.clear();
	APBinner<T>::_histogram.resize(Analyzer::posbins, Histogram_1D(Analyzer::angbins, 0));
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
		this->binner.HistogramOutput(output);
	  return;
	}

};
#endif
