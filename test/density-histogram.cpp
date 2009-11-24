#include "density-histogram.h"

template <class T>
void DensityBinner<T>::operator() (T * t) {
  // if procesing a new atom type that hasn't yet been encountered, create a new histogram for it
  if (!MapKeyExists(t))
	AddNewHistogram(t);

  int positionbin = Analyzer::PositionBin (t->Position()[y]);
  _histograms[t->Name()][positionbin]++;

  return;
}

template <class T>
// Output data from the histograms
void DensityBinner<T>::Output (FILE * output) {

  rewind (output);

  // first print out a header consisting of column headers
  fprintf (output, "%13s", "Position");	// The position column
  // column headers for each atom type (name)
  for (Histogram_it it = _histograms.begin(); it != _histograms.end(); it++) {
	fprintf (output, "%-13s", (*it).first.c_str());
  }
  fprintf(output, "\n");

  double position;
  // for every position in the system
  for (int pos = 0; pos < Analyzer::posbins; pos++) {

	// print the slab position
	position = double(pos)*Analyzer::posres + Analyzer::posmin;
	fprintf (output, "% 12.5f", position);

	// go through each atom type
	for (Histogram_it it = _histograms.begin(); it != _histograms.end(); it++) {
	  // and print the density value for the given position in it's own column
	  double density = double((*it).second[pos]) / Analyzer::timesteps;
	  fprintf (output, "% 12.5f", density);
	}
  fprintf(output, "\n");
  }
  return;
}

int main () {

  WaterSystemParams wsp ("atomic-densities.dat");

  DensityAnalyzer da (wsp);

  da.SystemAnalysis();

  return 0;
}
