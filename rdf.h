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
typedef std::map<NamePair, Histogram2D<double> > RDF2D_map;

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
    PrintParams ();
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

  void PrintParams () const {

    printf ("Atom pairs parsed from configuration file:\n");
    for (NamePairList::const_iterator it = name_pairs.begin(); it != name_pairs.end(); it++) {
      std::cout << it->first << " -- " << it->second << std::endl;
    }

    printf ("\n\nPerforming an RDF analysis using the following parameters:\n");
    printf ("\tRDF Parameters:\n\t\tMinimum: %f\n\t\tMaximum: %f\n\t\tResolution: %f\n",
	rdf_min, rdf_max, rdf_res);
    printf ("\tPosition Parameters:\n\t\tMinimum: %f\n\t\tMaximum: %f\n\t\tResolution: %f\n\n",
	position_min, position_max, position_res);

    return;
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





/**********************************************************************************************/ 
/*************** Abstract base class for processing RDFs in multiple dimensions ***************/
template <class T>
class RDFProcessor : public std::binary_function<T,T,bool> {
  public:
    RDFProcessor (const RDFParameters& params)
    {
      _params = params;
      _total_volume =  4.0/3.0 * M_PI * pow(_params.rdf_max, 3);
    }

    /* process the RDF and histogramming of the distance between the two atoms */
    void operator() (const T atom1, const T atom2);

    /* print the histograms to a file */
    virtual void Output (FILE * output) const = 0;

    static RDFParameters _params;
    static double _total_volume;	// the sphere of bonding distances in which we're calculating RDFs

    virtual void BinAtomPairData (const NamePair& name_pair, const T a1, const T a2) = 0;
};

template <class T> RDFParameters RDFProcessor<T>::_params;
template <class T> double RDFProcessor<T>::_total_volume;

  template <class T>
void RDFProcessor<T>::operator() (const T atom1, const T atom2)
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







/**********************************************************************************************/ 
/****************************** 1D RDF Processor **********************************************/ 
template <class T>
class RDFMachine : public RDFProcessor<T> {
  public:

    RDFMachine (const RDFParameters& params);

    virtual void Output (FILE * output) const;

  private:
    static RDF_map _rdfs;

    /* Calculates the rdf data and adds it into the correct histogram */
    virtual void BinAtomPairData (const NamePair& name_pair, const T a1, const T a2);

};

template <class T> RDF_map RDFMachine<T>::_rdfs;

  template <class T>
  RDFMachine<T>::RDFMachine (const RDFParameters& rdfparams)
: RDFProcessor<T> (rdfparams)
{

  // Create a new histogram for each pair of atomic names required for the analysis
  for (NamePairList::const_iterator it = RDFProcessor<T>::_params.begin(); it != RDFProcessor<T>::_params.end(); it++) {
    _rdfs.insert(std::make_pair(
	  *it,
	  Histogram1D<double>(RDFProcessor<T>::_params.rdf_min, RDFProcessor<T>::_params.rdf_max, RDFProcessor<T>::_params.rdf_res)));
  }

  printf ("RDF Analysis is using only one dimension - interatomic distances\n");
  return;
}

  template <class T> 
void RDFMachine<T>::BinAtomPairData (const NamePair& name_pair, const T a1, const T a2) 
{
  // The inter-atomic distance is calculated
  double atomic_distance = MDSystem::Distance(a1, a2).Magnitude();
  // and then the histogram is updated to track populations
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
  double min = RDFProcessor<T>::_params.rdf_min;
  double max = RDFProcessor<T>::_params.rdf_max;
  double res = RDFProcessor<T>::_params.rdf_res;
  int size = (int)((max - min)/res);
  double rdf_position;


  double dV, n, N;

  for (int rdf_i = 1; rdf_i < size; rdf_i++)
  {
    rdf_position = (double)rdf_i * res + min;
    // The differential volume in the spherical shell
    dV = 4.0 * M_PI * pow(rdf_position, 2) * res;

    // print out the rdf-position at the start of each row
    fprintf (output, "%13f", rdf_position);

    /* then output a data point for each rdf at that position */
    for (RDF_map::iterator it = _rdfs.begin(); it != _rdfs.end(); it++)
    {
      // the number of the given atom-pairs separated by a given distance
      n = double(it->second.Population(rdf_position));
      // Total number of pairs that were processed for the rdf
      N = double (it->second.Count());

      fprintf (output, "%13f", n * RDFProcessor<T>::_total_volume / dV / N);
    }
    fprintf (output, "\n");
  }

  return;
}


/************************* 2D RDF *****************************/
/* This functor takes atom-pairs to generate a radial-distribution-functions.
 * Application of the functor on an atom-pair creates several RDFs. Each atom
 * pair has its own set of RDFs, and each set of RDFs is divided up into sub-
 * RDFs representing slices of volume cutting up an MD slab.
 * The slab-position of an RDF is based on the first atom named in an atom-pair.
 * All the atom pairs are defined a priori and supplied to the functor's ctor. */
template <class T>
class RDF2DMachine : public RDFProcessor<T> {
  public:

    RDF2DMachine (const RDFParameters& params);

    /* Adds the given pair into the histogram after checking to see if that pair is one of those being analyzed */
    virtual void Output (FILE * output) const;

  private:
    static RDF2D_map _rdfs;

    /* Calculates the rdf data and adds it into the correct histogram */
    virtual void BinAtomPairData (const NamePair& name_pair, const T a1, const T a2);

};

template <class T> RDF2D_map RDF2DMachine<T>::_rdfs;

  template <class T>
  RDF2DMachine<T>::RDF2DMachine (const RDFParameters& params) 
: RDFProcessor<T>(params)
{
  // Create a new histogram for each pair of atomic names required for the analysis
  // ordering of all successive calls to the 2D histogram is [system position][rdf position/distance]
  for (NamePairList::const_iterator it = RDFProcessor<T>::_params.begin(); it != RDFProcessor<T>::_params.end(); it++) {
    _rdfs.insert(std::make_pair( 
	  *it, Histogram2D<double>( 
	    std::make_pair(RDFProcessor<T>::_params.position_min, RDFProcessor<T>::_params.rdf_min),
	    std::make_pair(RDFProcessor<T>::_params.position_max, RDFProcessor<T>::_params.rdf_max),
	    std::make_pair(RDFProcessor<T>::_params.position_res, RDFProcessor<T>::_params.rdf_res))));
  }

	printf ("RDF Analysis is using two dimensions - slab position & interatomic distances\n");

  return;
}

  template <class T> 
void RDF2DMachine<T>::BinAtomPairData (const NamePair& name_pair, const T a1, const T a2) 
{
  // The inter-atomic distance
  double atomic_distance = MDSystem::Distance(a1, a2).Magnitude();
  // Because the atomic position is being taken into account when slicing up the slab into many pieces, the correct atom must be used to determine the rdf position bin.
  // Thus, the following will choose the correct atomic position to use (namely, the first one listed in the name pair)
  double slab_position = (a1->Name() == name_pair.first) ? a1->Position()[WaterSystem<AmberSystem>::axis] : a2->Position()[WaterSystem<AmberSystem>::axis];
  // then call the 2Dhistogram to do the actual binning using the position and inter-atomic distances
  _rdfs.find(name_pair)->second(slab_position, atomic_distance);

  return;
}

template <class T>
void RDF2DMachine<T>::Output (FILE * output) const
{

  rewind(output);

  // first a header that lists Position, name-pairs, and min/max/res information for both position and rdf 
  // first column is the rdf position 
  fprintf (output, "%13s", "Position");
  // next the columns are for each name-pair
  for (RDF2D_map::iterator it = _rdfs.begin(); it != _rdfs.end(); it++)
  {
    std::string name1 = it->first.first;
    std::string name2 = it->first.second;
    fprintf (output, " (%s-%s)     ", name1.c_str(), name2.c_str());
  }
  fprintf (output, "\n");

  // to keep track of all the information of the system analysis, the next line will list 6 numbers:
  // position min, max, resolution, and rdf min, max, and resolution.
  // next the position extents/res 
  fprintf (output, "position %f %f %f ", RDFProcessor<T>::_params.position_min, RDFProcessor<T>::_params.position_max, RDFProcessor<T>::_params.position_res);
  // then the rdf information 
  fprintf (output, "rdf %f %f %f\n", RDFProcessor<T>::_params.rdf_min, RDFProcessor<T>::_params.rdf_max, RDFProcessor<T>::_params.rdf_res);

  // The differential volume in the spherical shell
  double dV, n, N;

  // main rdf-calculation loop to process populations into rdfs 
  // go through each rdf slice
  for (
      double rdf_distance = RDFProcessor<T>::_params.rdf_min; 
      rdf_distance < RDFProcessor<T>::_params.rdf_max; 
      rdf_distance += RDFProcessor<T>::_params.rdf_res)
  {
    // print out the rdf-position
    fprintf (output, "% 13f", rdf_distance);
    // calculate the differential shell volume
    dV = 4.0 * M_PI * pow(rdf_distance, 2) * RDFProcessor<T>::_params.rdf_res;

    // for each pair of atoms
    for (RDF2D_map::iterator it = _rdfs.begin(); it != _rdfs.end(); it++)
    {
      // for each slab position in the region of interest
      for (double pos = RDFProcessor<T>::_params.position_min; pos < RDFProcessor<T>::_params.position_max; pos += RDFProcessor<T>::_params.position_res)
      {
	// find the population of the given atom-pairs separated by a given distance
	n = double(it->second.Population(pos, rdf_distance));
	// Total number of pairs that were processed for the rdf
	N = double (it->second.Count(pos));

	// and print out the rdf-value for the given slab position, rdf-distance combo
	fprintf (output, " %8.4f", n * RDFProcessor<T>::_total_volume / dV / N);
      }
    }
    fprintf (output, "\n");
  }

  return;
}

#endif
