#include "rdf.h"

// given two atoms, the names and positions are checked, and then the histogram is updated 
void RDF::operator() (ap_t a, ap_t b) {
  double atomic_distance = MDSystem::Distance(a, b).Magnitude();
  if (CheckBounds(atomic_distance)) {
    UpdateHistogram(atomic_distance);
  }
  return;
}

void RDF::DataOutput (FILE * output) {
  rewind(output);

  // Header row
  fprintf (output, "Position %s %s", params.first.c_str(), params.second.c_str());
  fprintf (output, "\n");

  // Info about the extents for analyzing later and graphics output
  fprintf (output, "min:% 8.3f max:% 8.3f res:% 8.3f\n", params.min, params.max, params.res);

  /**** Data Rows ****/
  double population;

  printf ("%f\t%f\n", totalVolume, accessCount);
  for (double i = 0; i < histogram.size(); i++)
  {
    // the RDF position
    double pos = Position(i);
    // The differential volume in the spherical shell of radius pos
    shellVolume = 4.0 * M_PI * pow(pos, 2) * params.res;

    // print out the rdf-position at the start of each row
    fprintf (output, "% 13.3f", pos);

    // then output a data point for each position
    // the number of the given atom-pairs separated by a given distance = double(histogram[i]);
    // Total number of pairs that were processed for the rdf = double(accessCount);
    fprintf (output, "% 13.5f\n", histogram[i] * totalVolume / shellVolume / accessCount);
  }
  fflush(output);
  return;
}

void RDFAnalyzer::Setup () {
  LoadAll();
  std::pair<double,double> slice = std::make_pair (params.posmin, params.posmax);
  // slice the slab into the region of interest
  KeepAtomsInSlice(int_atoms, slice);

  // find all the unique atom names of atoms that we'll be processing through
  std::vector<std::string> atom_names;
  atom_names.push_back(params.first);
  atom_names.push_back(params.second);

  // keep only the atoms that appear in the list of names that we'll be processing
  KeepAtomsByNames(int_atoms, atom_names);
  KeepAtomsByNames(sys_atoms, atom_names);

  return;
}

void RDFAnalyzer::Analysis () {
  // Send each atom pair to the RDF machinery
  // perform the correlation between each atom in the specified interface/slice...
  for (Atom_it it = int_atoms.begin(); it != int_atoms.end(); it++)
  {
    // to every other atom in the system
    for (Atom_it jt = sys_atoms.begin(); jt != sys_atoms.end(); jt++)
    {
      if (*it == *jt) 
	continue;

      rdf(*it, *jt);
    }
  }
  return;
}

void RDFAnalyzer::DataOutput (const unsigned int timestep)
{ 
  rdf.DataOutput(output); 
  Setup();
  return;
}

int main () {

  libconfig::Config cfg;
  cfg.readFile("system.cfg");

  WaterSystemParams wsp (cfg);
  RDFAnalyzer rdf (wsp);
  rdf.SystemAnalysis();

  return 0;
}
