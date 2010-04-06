#ifndef __RDF_ANALYSIS_H
#define __RDF_ANALYSIS_H

#include "../utility.h"
#include "../analysis.h"
#include "../rdf.h"
#include <libconfig.h++>

// abstract class for RDF analyzers
template <class T>
class RDFAnalyzer : public Analyzer {

  public:
    RDFAnalyzer (const WaterSystemParams& wsp)
      : Analyzer(wsp), rdfparams(wsp.config_file), rdf(rdfparams)
    { return; }

  protected:
    
    virtual void Setup ();
    virtual void Analysis ();
    virtual void PostAnalysis ();
    virtual void DataOutput (const unsigned int timestep);

    // the utility functor for getting all the data accumulated
    RDFParameters rdfparams;
    T rdf;
};

#endif
