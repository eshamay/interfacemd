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
typedef std::map<NamePair_t, Histogram2D<double> > RDF_map;
	
/* This functor takes atom-pairs to generate a radial-distribution-functions.
 * Application of the functor on an atom-pair creates several RDFs. Each atom
 * pair has its own set of RDFs, and each set of RDFs is divided up into sub-
 * RDFs representing slices of volume cutting up an MD slab.
 * The slab-position of an RDF is based on the first atom named in an atom-pair.
 * All the atom pairs are defined a priori and supplied to the functor's ctor. */
template <class T>
class RDFMachine : public std::binary_function<T,T,bool> {
  public:
	// Initialization of all the histograms - the parameter pairs refer to <slab position, rdf>
	RDFMachine (const NamePairList& names, const double_pair maxima, const double_pair minima, const double_pair bin_widths);

	/* Adds the given pair into the histogram after checking to see if that pair is one of those being analyzed */
	void operator() (const T a1, const T a2) 
	{
	  /* Check if the pair's RDF is being calculated i.e. in the list of name-pairs */
	  NamePair_t name_pair = std::make_pair(a1->Name(), a2->Name());
	  /* find the actual way the pair is ordered in the list (i.e. first atom named = first atom, etc) */
	  NamePairList::iterator list_pair = PairListMember(name_pair, _name_pair_list.begin(), _name_pair_list.end());

	  if (list_pair != _name_pair_list.end())	// see if the atomic pair is found in the list
		BinAtomPairData (*list_pair, a1, a2);

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
RDFMachine<T>::RDFMachine (const NamePairList& names, const double_pair maxima, const double_pair minima, const double_pair bin_widths)
{
  // Create a new histogram for each pair of atomic names required for the analysis
  RUN (names) {
	//_rdfs[names[i]] = Histogram2D<double>(maxima, bin_widths, minima);
	_rdfs.insert(std::make_pair(names[i], Histogram2D<double>(maxima, bin_widths, minima)));
	_name_pair_list.push_back (names[i]);
  }
  return;
}

template <class T> 
void RDFMachine<T>::BinAtomPairData (const NamePair_t name_pair, const T a1, const T a2) 
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

/* Print the data from the RDFs to the given output file */
template <class T> 
void RDFMachine<T>::Output (FILE * output) const
{
  /* first print out a header row that contains the position and pair-names */
  rewind(output);

  std::pair<int,int> size = make_pair(0,0);
  std::pair<double,double> resolution = make_pair(0.0,0.0);
  std::pair<double,double> max = make_pair(0.0,0.0);
  std::pair<double,double> min = make_pair(0.0,0.0);

  fprintf (output, "%13s", "Distance");
  for (RDF_map::iterator it = _rdfs.begin(); it != _rdfs.end(); it++)
  {
	if (!size.first) 
	{
	  size = it->second.size;
	  resolution = it->second.resolution;
	  max = it->second.max;
	  min = it->second.min;
	}

	std::string name1 = it->first.first;
	std::string name2 = it->first.second;
	for (double slab_position = min.first; slab_position < max.first; slab_position += resolution.first)
	{
	  // output the atom names in the pair
	  fprintf (output, "(%2s-%2s)%4.2f ", name1.c_str(), name2.c_str(), slab_position);
	}
  }
  fprintf (output, "\n");

  /* print out on each row the rdf position */
  for (int rdf_i = 0; rdf_i < size.second; rdf_i++)
  {
	double rdf_position = double(rdf_i)*resolution.second;	// assuming that all the RDFs will have the same resolution

	// print out the rdf-position at the start of each row
	fprintf (output, "%13f", rdf_position);		// don't add the rdf min because it will always be 0... for now.

	/* then output a data point for each atom name-pair's rdf value for each slab position */
	for (RDF_map::iterator it = _rdfs.begin(); it != _rdfs.end(); it++)
	{
	  size = it->second.size;
	  resolution = it->second.resolution;
	  min = it->second.min;

	  /* and for each name-pair, output the columns for each of the slab positions */
	  for (int slab_i = 0; slab_i < size.first; slab_i++)
	  {
		// the number density of the given atom-pair separated by a given distance
		double num = double(it->second.Element(slab_i, rdf_i));
		// the differential volume of the shell being parsed
		double volume = 4.0/3.0*M_PI*(pow(rdf_position, 3));
		// here's our normalization
		double rdf = volume / (4 * M_PI * pow(rdf_position, 2) * resolution.second) / num;	

		fprintf (output, "%13f", rdf);
	  }
	}
	fprintf (output, "\n");
  }

  return;
}


#endif
