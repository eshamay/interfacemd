#ifndef __RDF_ANALYSIS_H
#define __RDF_ANALYSIS_H

#include "../utility.h"
#include "../analysis.h"
#include "../rdf.h"

class RDFAnalyzer : public Analyzer {

  public:
	RDFAnalyzer (WaterSystemParams& wsp, const NamePairList& names, const double_pair maxima, const double_pair minima, const double_pair bin_widths) 
	  : Analyzer(wsp), rdf(names, maxima, minima, bin_widths)
	{ return; }

	// the utility functor for getting all the data accumulated
	RDFMachine<Atom *> rdf;

	void Setup () 
	{ 
	  LoadAll();
	  SLICE_BY_POSITION(int_atoms, Atom *, 28.0, 29.0);
	  return;
	}

	void Analysis () {
	  /* Send each atom pair to the RDF machinery */
	  for (int i = 0; i < int(int_atoms.size()) - 1; i++)
	  {
		for (int j = i+1; j < int(int_atoms.size()); j++)
		{
		  rdf(int_atoms[i], int_atoms[j]);
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
