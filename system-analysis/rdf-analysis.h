#ifndef __RDF_ANALYSIS_H
#define __RDF_ANALYSIS_H

#include "../utility.h"
#include "../analysis.h"
#include "../rdf.h"
#include <libconfig.h++>

class RDFAnalyzer : public Analyzer {

  public:
    RDFAnalyzer (const WaterSystemParams& wsp)
      : Analyzer(wsp), rdfparams(wsp.config_file), rdf(rdfparams)
    { return; }


  private:
    // the utility functor for getting all the data accumulated
    RDFParameters rdfparams;
    RDFMachine<Atom *> rdf;
    NamePairList _atom_name_pairs;	// names of atom pairs that will be processed

    void Setup ();
    void Analysis ();
    void PostAnalysis ();
    void DataOutput (const unsigned int timestep);

};

#endif
