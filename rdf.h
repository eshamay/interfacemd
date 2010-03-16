#ifndef __RDF_H_
#define __RDF_H_

#include <map>
#include <utility>
#include "utility.h"
#include "atom.h"

typedef std::vector<int> Histogram_t;
typedef std::pair<std::string, std::string> NamePair;
typedef std::vector<NamePair> NamePairList;
typedef std::pair<double, double> double_pair;
typedef std::map<NamePair, Histogram1D<double> > RDF_map;

struct RDFParameters {

  RDFParameters () { return; }

  /* Parse the RDF parameters from the given configuration file */
  RDFParameters (libconfig::Config * config_file)
  {
    try {
      rdf_min = config_file->lookup("analysis.rdf.minimum");
      rdf_max = config_file->lookup("analysis.rdf.maximum");
      rdf_res = config_file->lookup("analysis.rdf.resolution");
      position_min = config_file->lookup("analysis.rdf.position-cutoff-low");
      position_max = config_file->lookup("analysis.rdf.position-cutoff-high");
      position_res = config_file->lookup("analysis.rdf.position-resolution");
    }
    catch(const libconfig::SettingTypeException &stex) {
      std::cerr << "Something is wrong with the configuration parameters or file - check syntax\n(watersystem.h)" << std::endl;
      exit(EXIT_FAILURE);
    }
    catch(const libconfig::SettingNotFoundException &snfex) {
      std::cerr << "Missing one of the rdf settings in the configuration file" << std::endl;
      exit(EXIT_FAILURE);
    }

    ParseAtomPairs (config_file);
  }

  // Parse the atom pairs from the configuration file
  void ParseAtomPairs (libconfig::Config * cfg) {

    try
    {
      libconfig::Setting &atompairs = cfg->lookup("analysis.rdf.atom-pairs");
      for (int i = 0; i < atompairs.getLength(); i++)
      {
	std::string atom1 = atompairs[i][0];
	std::string atom2 = atompairs[i][1];
	name_pairs.push_back(std::make_pair(atom1,atom2));
      }
    }
    catch(const libconfig::SettingTypeException &stex) {
      std::cerr << "Check the atom-pairs setting in the RDF configuration section" << std::endl;
      exit(EXIT_FAILURE);
    }
    catch(const libconfig::SettingNotFoundException &snfex) {
      std::cerr << "Missing the atom-pairs setting (analysis.rdf.atom-pairs) in the configuration file" << std::endl;
      exit(EXIT_FAILURE);
    }
    return;
  }

  void PrintParams () {

    printf ("Atom pairs parsed from configuration file:\n");
    for (NamePairList::iterator it = name_pairs.begin(); it != name_pairs.end(); it++) {
      std::cout << it->first << " -- " << it->second << std::endl;
    }

    printf ("rdf_max = %f\nrdf_min = %f\nrdf_res = %f\n\n", rdf_max, rdf_min, rdf_res);
    printf ("position_max = %f\nposition_min = %f\nposition_res = %f\n\n", position_max, position_min, position_res);
  }

  typedef NamePairList::const_iterator namepair_it;
  namepair_it begin () { return name_pairs.begin(); }
  namepair_it end () { return name_pairs.end(); }

  double rdf_min;
  double rdf_max;
  double rdf_res;
  double position_min;
  double position_max;
  double position_res;
  NamePairList name_pairs;

};


/* This functor takes atom-pairs to generate a radial-distribution-functions.
 * Application of the functor on an atom-pair creates several RDFs. Each atom
 * pair has its own set of RDFs, and each set of RDFs is divided up into sub-
 * RDFs representing slices of volume cutting up an MD slab.
 * The slab-position of an RDF is based on the first atom named in an atom-pair.
 * All the atom pairs are defined a priori and supplied to the functor's ctor. */
template <class T>
class RDFMachine : public std::binary_function<T,T,bool> {
  public:

    RDFMachine (const RDFParameters& params);

    /* Adds the given pair into the histogram after checking to see if that pair is one of those being analyzed */
    void operator() (const T atom1, const T atom2) 
    {

      /** Check if the pair's RDF is being calculated i.e. in the list of name-pairs **/
      NamePair name_pair = std::make_pair(atom1->Name(), atom2->Name());

      /* find the actual way the pair is ordered in the list (i.e. first atom named = first atom, etc) */
      NamePairList::iterator list_pair = PairListMember(
	  name_pair, _params.name_pairs.begin(), _params.name_pairs.end());

      if (list_pair < _params.name_pairs.end())	// see if the atomic pair is found in the list
      {
	BinAtomPairData (*list_pair, atom1, atom2);
      }
      return;
    }

    /* print the histograms to a file */
    void Output (FILE * output) const;

  private:
    static RDFParameters _params;
    static unsigned long _total_atoms;
    static RDF_map _rdfs;
    static double _total_volume;	// the sphere of bonding distances in which we're calculating RDFs

    /* Calculates the rdf data and adds it into the correct histogram */
    void BinAtomPairData (const NamePair& name_pair, const T a1, const T a2);

};

template <class T> RDFParameters RDFMachine<T>::_params;
template <class T> RDF_map RDFMachine<T>::_rdfs;
template <class T> double RDFMachine<T>::_total_volume;

  template <class T>
RDFMachine<T>::RDFMachine (const RDFParameters& rdfparams)
{
  _params = rdfparams;

  _total_volume =  4.0/3.0 * M_PI * pow(_params.rdf_max, 3);

  // Create a new histogram for each pair of atomic names required for the analysis
  for (NamePairList::const_iterator it = _params.begin(); it != _params.end(); it++) {
    //RUN (_params.name_pairs) {
    _rdfs.insert(std::make_pair(
	  *it,
	  Histogram1D<double>(_params.rdf_min, _params.rdf_max, _params.rdf_res)));
  }

  return;
}


  template <class T> 
void RDFMachine<T>::BinAtomPairData (const NamePair& name_pair, const T a1, const T a2) 
{
  // The inter-atomic distance
  double atomic_distance = MDSystem::Distance(a1, a2).Magnitude();
  // Because the atomic position is being taken into account when slicing up the slab into many pieces, the correct atom must be used to determine the rdf position bin.
  // Thus, the following will choose the correct atomic position to use (namely, the first one listed in the name pair)
  //double slab_position = (a1->Name() == name_pair.first) ? a1->Position()[WaterSystem<AmberSystem>::axis] : a2->Position()[WaterSystem<AmberSystem>::axis];
  // then call the 2Dhistogram to do the actual binning using the position and inter-atomic distances
  //_rdfs.find(name_pair)->second(slab_position, atomic_distance);
  _rdfs.find(name_pair)->second(atomic_distance);
  return;
}

/* Print the data from the RDFs to the given output file */
template <class T> 
void RDFMachine<T>::Output (FILE * output) const
{
  rewind(output);

  /* first print out a **** header row **** that contains the position and pair-names */
  /* first column is the rdf position */
  fprintf (output, "%13s", "Position");
  /* run through each histogram (atom-pair) and print out the header information */
  for (RDF_map::iterator it = _rdfs.begin(); it != _rdfs.end(); it++)
  {
    std::string name1 = it->first.first;
    std::string name2 = it->first.second;
    fprintf (output, " (%s-%s)     ", name1.c_str(), name2.c_str());
  }
  fprintf (output, "\n");

  /**** Data Rows ****/
  /* print out on each row the rdf position - assuming all the RDFs have the same min,max,res... */
  double min = _params.rdf_min;
  double max = _params.rdf_max;
  double res = _params.rdf_res;
  int size = (int)((max - min)/res);
  double rdf_position;


  for (int rdf_i = 1; rdf_i < size; rdf_i++)
  {
    rdf_position = (double)rdf_i * res + min;

    // print out the rdf-position at the start of each row
    fprintf (output, "%13f", rdf_position);

    /* then output a data point for each rdf at that position */
    for (RDF_map::iterator it = _rdfs.begin(); it != _rdfs.end(); it++)
    {
      // The differential volume in the spherical shell
      double dV = 4.0 * M_PI * pow(rdf_position, 2) * res;
      // the number of the given atom-pairs separated by a given distance
      double n = double(it->second.Population(rdf_position));
      // Total number of pairs that were processed for the rdf
      double N = double (it->second.Count());

      fprintf (output, "%13f", n * _total_volume / dV / N);
    }
    fprintf (output, "\n");
  }

  return;
}

#endif
