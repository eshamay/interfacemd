#ifndef __RDF_ANALYSIS_H
#define __RDF_ANALYSIS_H

#include "../utility.h"
#include "../analysis.h"
#include "../rdf.h"

class RDFAnalyzer : public Analyzer {

  public:
	RDFAnalyzer (WaterSystemParams& wsp, const NamePairList& names, const double_pair maxima, const double_pair minima, const double_pair bin_widths) 
	  : Analyzer(wsp), rdf(names, maxima, minima, bin_widths), _maxima(maxima), _minima(minima), _atom_name_pairs(names)
	{ return; }

	// the utility functor for getting all the data accumulated
	RDFMachine<Atom *> rdf;
	double_pair _maxima;
	double_pair _minima;
	NamePairList _atom_name_pairs;

	void Setup () 
	{ 
	  LoadAll();
	  // slice the slab into the region of interest
	  SLICE_BY_POSITION(int_atoms, Atom *, _minima.first, _maxima.first);

	  // find all the unique atom names of atoms that we'll be processing through
	  std::vector<string> atom_names;
	  UNIQUE_PAIRLIST_ELEMENTS(_atom_name_pairs, atom_names);
	  // keep only the atoms that appear in the list of names that we'll be processing
	  KEEP_BY_NAMES(int_atoms, Atom *, atom_names);

	  return;
	}

	void Analysis () {
	  /* Send each atom pair to the RDF machinery */
	  for (Atom_ptr_vec::iterator it = int_atoms.begin(); it < int_atoms.end(); it++)
	  {
		for (Atom_ptr_vec::iterator jt = int_atoms.begin(); jt < int_atoms.end(); jt++)
		{
		  rdf(*it, *jt);
		}
	  }
	  return;
	}

	// nothing at the moment
	void PostAnalysis () { return; }

	void DataOutput (const unsigned int timestep) {
	  if (!(timestep % (output_freq * 10)))
		rdf.Output(output);
	  return;
	}
};

#endif
