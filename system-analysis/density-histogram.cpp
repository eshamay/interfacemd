#include "density-histogram.h"

template <class T, class U>
void DensityBinner<T,U>::operator() (T t) {
  // if procesing a new atom type that hasn't yet been encountered, create a new histogram for it
  if (!MapKeyExists(t))
	AddNewHistogram(t);

  _histograms[t->Name()][Analyzer<U>::PositionBin(t)]++;

  return;
}

template <class T, class U>
// Output data from the histograms
void DensityBinner<T,U>::Output (FILE * output, const int timestep) {

  rewind (output);

  double scale = (!timestep) ? 1.0 : double(timestep);

  // first print out a header consisting of column headers
  fprintf (output, "Position");	// The position column
  // column headers for each atom type (name)
  for (Histogram_it it = _histograms.begin(); it != _histograms.end(); it++) {
	fprintf (output, "%13s", (*it).first.c_str());
  }
  fprintf(output, "\n");

  double position;
  // for every position in the system
  for (int pos = 0; pos < Analyzer<U>::posbins; pos++) {

	// print the slab position
	position = double(pos)*Analyzer<U>::posres + Analyzer<U>::posmin;
	fprintf (output, "% 13.5f", position);

	// go through each atom type
	for (Histogram_it it = _histograms.begin(); it != _histograms.end(); it++) {
	  // and print the density value for the given position in it's own column
	  double density = double((*it).second[pos]) / scale;
	  fprintf (output, "% 13.5f", density);
	}
  fprintf(output, "\n");
  }
  return;
}

int main () {

  libconfig::Config cfg;
  cfg.readFile("system.cfg");

  std::string filename = cfg.lookup("analysis.density.filename");
  libconfig::Setting &analysis = cfg.lookup("analysis");
  analysis.add("filename", libconfig::Setting::TypeString) = filename;

  WaterSystemParams wsp (cfg);

  DensityAnalyzer<GMXSystem> da (wsp);

  da.SystemAnalysis();

  return 0;
}
