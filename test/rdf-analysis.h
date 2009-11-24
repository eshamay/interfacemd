#ifndef __RDF_ANALYSIS_H
#define __RDF_ANALYSIS_H

#include "../utility.h"
#include "../analysis.h"
#include "../rdf.h"

class RDFAnalyzer : public Analyzer {

  public:
	RDFAnalyzer (WaterSystemParams& wsp, const NamePairList& names, const double max_length, const double bin_size) 
	  : Analyzer(wsp), rdf(names, max_length, bin_size) 
	{ return; }

	// the utility functor for getting all the data accumulated
	RDFMachine<Atom *> rdf;

	void Setup () 
	{ 
	  FindAll (); 
	  return;
	}

	void Analysis () {
	  /* Send each atom pair to the RDF machinery */
	  for (int i = 0; i < int_atoms.size() - 1; i++)
	  {
		for (int j = i+1; j < int_atoms.size(); j++)
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
