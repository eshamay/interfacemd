#include "efield.h"

template <class T>
void EFieldBinner<T>::operator() (Atom * t) {

	charge = AtomCharge(t);
	position = t->Position()[WaterSystem<T>::axis];

	// for each atom being analyzed, go through each position from +/- 10.0 angstroms (potential cutoff) and bin the contribution to the Efield that it gives
	for (double z = position - cutoff; z < position + cutoff; z += Analyzer<T>::posres) {

		// only calculate the field within the slab-slice of interest
		if (z < WaterSystem<T>::posmin || z > WaterSystem<T>::posmax) continue;
		
		// find the distance from the atom to the point of calculation
		distance = (position > z) ? position-z : z-position;

		if (distance < 0.1) continue;

		// calculate the electric field component in the major-axis direction
		efield = charge/(distance*distance);
		// bin the electric field contribution by summing it in the histogram for the system position
		_histogram[int((z-WaterSystem<T>::posmin)/Analyzer<T>::posres)] += efield;
	}

  return;
}

template <class U>
double EFieldBinner<U>::AtomCharge (Atom * t) {
  double ret;

  if (t->Name().find("H") != std::string::npos) ret = 0.365;
  //else if (t->Name().find("O1") != std::string::npos) ret = -0.65;
  //else if (t->Name().find("O2") != std::string::npos) ret = -0.65;
  //else if (t->Name().find("O3") != std::string::npos) ret = -0.65;
  else if (t->Name().find("S") != std::string::npos) ret = 2.0;
  else if (t->Name().find("O1") != std::string::npos) ret = -1.0;
  else if (t->Name().find("O2") != std::string::npos) ret = -1.0;
  else if (t->Name().find("O3") != std::string::npos) ret = -1.0;
  else if (t->Name().find("O4") != std::string::npos) ret = -1.0;
  else if (t->Name().find("O") != std::string::npos) ret = -0.73;
  else if (t->Name().find("Cl1") != std::string::npos) ret = 0.0404;
  else if (t->Name().find("Cl2") != std::string::npos) ret = 0.0404;
  else if (t->Name().find("Cl3") != std::string::npos) ret = 0.0404;
  else if (t->Name().find("Cl4") != std::string::npos) ret = 0.0404;
  else if (t->Name().find("Cl") != std::string::npos) ret = -1.0;
  else if (t->Name().find("C") != std::string::npos) ret = -0.1616;
  else if (t->Name().find("NA") != std::string::npos) ret = 1.0;
  else if (t->Name().find("N") != std::string::npos) ret = 0.95;

  return ret;
}

template <class U>
// Output data from the histograms
void EFieldBinner<U>::Output (FILE * output, const int timestep) {

  rewind (output);

  double scale = (!timestep) ? 1.0 : double(timestep);

  // first print out a header consisting of column headers
  fprintf (output, "Position\tEField\n");	// The position column

  double position;
  double field;
  // for every position in the system
  for (int pos = 0; pos < Analyzer<U>::posbins; pos++) {
    // print the slab position
    position = double(pos)*Analyzer<U>::posres + Analyzer<U>::posmin;
    field = _histogram[pos] / scale;
    fprintf (output, "% 13.5f% 13.5f\n", position, field);
  }
return;
}

int main () {

  libconfig::Config cfg;
  cfg.readFile("system.cfg");

  std::string filename = cfg.lookup("analysis.efield.filename");
  libconfig::Setting &analysis = cfg.lookup("analysis");
  analysis.add("filename", libconfig::Setting::TypeString) = filename;

  WaterSystemParams wsp (cfg);

  EFieldAnalyzer<AmberSystem> ea (wsp);

  ea.SystemAnalysis();

  return 0;
}
