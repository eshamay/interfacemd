#include "rdf-analysis.h"

int main () {

  WaterSystemParams wsp ("system-RDFs.dat", 1);

  NamePairList names;
  names.push_back(std::make_pair("O", "O"));
  names.push_back(std::make_pair("O", "H1"));
  names.push_back(std::make_pair("O", "H2"));
  names.push_back(std::make_pair("O", "NA"));
  names.push_back(std::make_pair("O", "N"));

  RDFAnalyzer rdf (wsp, names, make_pair(50.0,10.0), make_pair(20.0,0.0), make_pair(2.0,0.025));

  rdf.SystemAnalysis();

  return 0;
}
