#ifndef __RDF_H_
#define __RDF_H_

#include <map>
#include <utility>
#include "utility.h"
#include "atom.h"

typedef std::vector<int> Histogram_t;
typedef std::pair<std::string, std::string> NamePair_t;
typedef std::pair<double, double> double_pair;
typedef std::vector<NamePair_t> NamePairList;
typedef std::map<NamePair_t, 2DHistogram> RDF_map;
	
/* This functor takes atom-pairs to generate a radial-distribution-functions.
 * Application of the functor on an atom-pair creates several RDFs. Each atom
 * pair has its own set of RDFs, and each set of RDFs is divided up into sub-
 * RDFs representing slices of volume cutting up an MD slab.
 * The slab-position of an RDF is based on the first atom named in an atom-pair.
 * All the atom pairs are defined a priori and supplied to the functors ctor. */
template <class T>
class RDFMachine : public std::binary_function<T,T,bool> {
  public:
	// the first values refer to the slab positions, and the 2nd value refer to the rdf itself
	RDFMachine (const NamePairList& names, double_pair maxima, double_pair minima, double_pair bin_widths);

	/* Adds the given pair into the histogram after checking to see if that pair is one of those being analyzed */
	void operator() (const T a1, const T a2) 
	{
	  /* Check if the pair's RDF is being calculated i.e. in the list of name-pairs */
	  NamePair_t name_pair = std::make_pair(a1->Name(), a2->Name());
	  /* find the actual way the pair is ordered in the list (i.e. first atom named = first atom, etc) */
	  NamePairList::iterator list_pair = PairListMember(name_pair, _name_pair_list.begin(), _name_pair_list.end());

	  if (list_pair != _name_pair_list.end())	// see if the atomic pair is found in the list
	  {
		BinAtomPairData (*list_pair, a1, a2);
	  }
	  return;
	}

	/* print the histograms to a file */
	void Output (FILE * output) const;

  private:
	static RDF_map _rdfs;
	static vector<NamePair_t> _name_pair_list;

	/* Calculates the rdf data and adds it into the correct histogram */
	void BinAtomPairData (const NamePair_t name_pair, const T a1, const T a2);

};

template <class T> RDF_map RDFMachine<T>::_rdfs;
template <class T> NamePairList RDFMachine<T>::_name_pair_list;

template <class T> 
RDFMachine<T>::RDFMachine (const NamePairList& names, double_pair maxima, double_pair minima, double_pair bin_widths)
{
  RDFMachine::_max_rdf_distance = max_length;
  RDFMachine::_bin_size = bin_size;

  // Create a new histogram for each pair of names for the analysis
  RUN (names) {
	_rdfs[names[i]] = 2DHistogram(maxima, bin_widths, minima);
	_name_pair_list.push_back (names[i]);
  }
  return;
}

template <class T> 
void RDFMachine<T>::BinAtomPairData (const NamePair_t name_pair, const T a1, const T a2) 
{
  // The inter-atomic distance
  double distance = MDSystem::Distance(a1, a2).Magnitude();
  // Because the atomic position is being taken into account when slicing up the slab into many pieces, the correct atom must be used to determine the rdf position bin.
  // Thus, the following will choose the correct atomic position to use (namely, the first one listed in the name pair)
  double position = (a1->Name() == name_pair.first) ? a1->Position()[WaterSystem::axis] : a2->Position()[WaterSystem::axis];
  // then call the 2Dhistogram to do the actual binning using the position and inter-atomic distances
  _rdfs[name_pair](position, distance);
  return;
}

/* Print the data from the RDFs to the given output file */
template <class T> 
void RDFMachine<T>::Output (FILE * output) const
{
  /* first print out a header row that contains the position and pair-names */
  rewind(output);

  fprintf (output, "Distance");
  for (RDF_map::iterator it = _rdfs.begin(); it != _rdfs.end(); it++)
  {
	fprintf (output, " (%-4s,%-4s) ", (*it).first.first.c_str(), (*it).first.second.c_str());	// output the atom names in the pair
  }
  fprintf (output, "\n");

  /* then go through each histogram and print out the data for the given position/distance */
  for (int i = 0; i < int(_max_rdf_distance/_bin_size + 1); i++)
  {
	fprintf (output, "%13f", double(i) * _bin_size);

	for (RDF_map::iterator it = _rdfs.begin(); it != _rdfs.end(); it++)
	{
	  double num = double(it->second[i]);					// the number density of the given atom-pair at a given distance
	  double volume = 4.0/3.0*M_PI*pow(_max_rdf_distance,3);		// total volume over which we're checking
	  double norm = volume / (4 * M_PI * pow(double(i)*_bin_size,2) * _bin_size) / num;	// here's our normalization

	  fprintf (output, "%13f", norm);
	}
	fprintf (output, "\n");
  }
  return;
}


#endif
