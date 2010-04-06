#include "rdf-analysis.h"

template <class T>
void RDFAnalyzer<T>::Setup () {
  LoadAll();
  // slice the slab into the region of interest
  SLICE_BY_POSITION(int_atoms, Atom *, rdfparams.position_min, rdfparams.position_max);

  // find all the unique atom names of atoms that we'll be processing through
  std::vector<string> atom_names;
  UNIQUE_PAIRLIST_ELEMENTS(rdfparams.name_pairs, atom_names);

  // keep only the atoms that appear in the list of names that we'll be processing
  KEEP_BY_NAMES(int_atoms, Atom *, atom_names);

  return;
}

template <class T>
void RDFAnalyzer<T>::Analysis () {
  /* Send each atom pair to the RDF machinery */
  for (Atom_it it = int_atoms.begin(); it != int_atoms.end(); it++)
  {
    for (Atom_it jt = it; jt != int_atoms.end(); jt++)
    {
      if (it == jt) 
	continue;

      rdf(*it, *jt);
    }
  }
  return;
}

template <class T>
void RDFAnalyzer<T>::PostAnalysis () { 
  return; 
}

template <class T>
void RDFAnalyzer<T>::DataOutput (const unsigned int timestep) {
  if (!(timestep % (output_freq * 10)))
  {
    rdf.Output(output);
    this->Setup();
  }
  return;
}

int main () {

  libconfig::Config cfg;
  cfg.readFile("system.cfg");

  /* 1d or 2d rdf */
  int rdf_type = cfg.lookup("analysis.rdf.rdf-type");

  /* get the output filename */
  std::string filename = cfg.lookup("analysis.rdf.filename");
  libconfig::Setting &analysis = cfg.lookup("analysis");
  analysis.add("filename", libconfig::Setting::TypeString) = filename;

  WaterSystemParams wsp (cfg);

  if (rdf_type == 1)
  {
  	RDFAnalyzer< RDFMachine<Atom *> > rdf (wsp);
  	rdf.SystemAnalysis();
  }
  else if (rdf_type == 2)
  {
  	RDFAnalyzer< RDF2DMachine<Atom *> > rdf (wsp);
  	rdf.SystemAnalysis();
  }

  return 0;
}
