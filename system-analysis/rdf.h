#ifndef __RDF_H
#define __RDF_H

#include <string>
#include "../analysis.h"
#include "../utility.h"

class RDFParameters {
  public:
    RDFParameters (const double _min, const double _max, const double _res, 
	const std::string _first, const std::string _second,
	const double _posmin, const double _posmax, const double _posres)
    : 
      min(_min), max(_max), res(_res), 
      first(_first), second(_second),
      posmin(_posmin), posmax(_posmax), posres(_posres)
    { return; }

    double min, max, res;
    double posmin, posmax, posres;
    std::string first, second;

    bool NamesMatch (const std::string a, const std::string b) const {
      return ((a == first && b == second) || (a == second && b == first));
    }

    bool InBounds (const double pos) const { 
      return (pos > min && pos < max);
    }
};
      

class RDF {
  private:

    typedef double storage_t;
    typedef std::vector<storage_t> histogram_t;
    typedef const Atom * ap_t;

    RDFParameters params;
    double totalVolume;
    histogram_t histogram;
    storage_t accessCount;
    double shellVolume;

    bool CheckNames (ap_t a, ap_t b) const { return params.NamesMatch (a->Name(), b->Name()); }
    bool CheckBounds (const double pos) const { return params.InBounds(pos); }

    // given a position, return the bin
    int Bin (const double pos) const {
      return int((pos - params.min)/params.res);
    }
    // given a bin, return the position
    double Position (const int bin) const {
      return double(bin)*params.res + params.min;
    }

    void UpdateHistogram (const double pos) {
      histogram[Bin(pos)]++;
      accessCount++;
    }

  public:
    RDF (const RDFParameters& pars) 
      : 
	params(pars),
	totalVolume (4.0/3.0 * M_PI * pow(params.max, 3)),
	// set up the histogram size/number of bins, and initialize them to 0
	histogram(int(params.max-params.min)/params.res, 0),
	accessCount(0)
    { return; }

    void operator() (ap_t a, ap_t b);
    void DataOutput (FILE * output);
};

class RDFAnalyzer : public Analyzer {
  public:
    RDFAnalyzer (const WaterSystemParams& wsp)
      : 
	Analyzer(wsp),
	params
	( wsp.config_file->lookup("analysis.rdf.minimum"),
	  wsp.config_file->lookup("analysis.rdf.maximum"),
	  wsp.config_file->lookup("analysis.rdf.resolution"),
	  wsp.config_file->lookup("analysis.rdf.atom-pairs")[0][0],
	  wsp.config_file->lookup("analysis.rdf.atom-pairs")[0][1],
	  wsp.config_file->lookup("analysis.rdf.position-cutoff-low"),
	  wsp.config_file->lookup("analysis.rdf.position-cutoff-high"),
	  wsp.config_file->lookup("analysis.rdf.position-resolution")),
	rdf(params)
  { return; }

  protected:

    RDFParameters params;
    RDF rdf;

    virtual void Setup ();
    virtual void Analysis ();
    virtual void DataOutput (const unsigned int timestep);

};

#endif
