#include "rdf-analysis.h"

int main () {

  WaterSystemParams wsp ("system-RDFs.dat");

  NamePairList names;
  names.push_back(std::make_pair("O", "O"));
  /*
  names.push_back(std::make_pair("O", "H"));
  names.push_back(std::make_pair("H", "H"));
  names.push_back(std::make_pair("O", "NA"));
  names.push_back(std::make_pair("H", "NA"));
  names.push_back(std::make_pair("N", "NA"));
  names.push_back(std::make_pair("N", "O"));
  names.push_back(std::make_pair("N", "H"));
  names.push_back(std::make_pair("C", "C"));
  */

  RDFAnalyzer rdf (wsp, names, make_pair(30.0,15.0), make_pair(29.0,0.0), make_pair(0.1,0.1));

  rdf.SystemAnalysis();

  return 0;
}
