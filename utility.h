#ifndef __UTIL_H
#define __UTIL_H

#include <vector>
#include <algorithm>
#include <utility>
#include <functional>
#include <iterator>
#include <iostream>
#include <cctype>
#include <cstdio>

namespace md_utility {

  /* A functor that takes a std::pair as a constructor argument, and all applications of the tester will test pairs against the initial one given. */
  template <class T>
	class EqualPairs : public std::binary_function<T,T,bool> {
	  public:
		bool operator() (const T& lhs, const T& rhs) const 
		{ return lhs == rhs || lhs.first == rhs.second && lhs.second == rhs.first; }
	};



  /* Searches for the supplied pair p in the sequence of pairs from first to last */
  template <class Iter> 
	Iter PairListMember (const typename std::iterator_traits<Iter>::value_type& p, Iter first, Iter last)
	{
	  return find_if (first, last, std::bind1st(EqualPairs<typename std::iterator_traits<Iter>::value_type>(), p));
	}


  // used for sorting/comparing pairs based on the first element
  template <typename T>
	class pair_sort_first_pred : public std::binary_function<T,T,bool> {
	  public:
		bool operator()(const T& left, const T& right) const {
		  return left.first < right.first;
		}
	};



  // sorts a container of pairs by the first element of the pairs
  template <typename Iter>
	void pair_sort_first (Iter first, Iter last) {
	  typedef typename std::iterator_traits<Iter>::value_type val_t;
	  std::sort (first, last, pair_sort_first_pred<val_t>());
	}	// Pair sort first





  // predicate to test if the name of an atom or molecule (determined by the template parameter) is found in a vector of names
  template <class U>
	class NameInList : public std::binary_function<U, std::vector<std::string>,bool> {
	  public:
		bool operator() (const U u, const std::vector<std::string>& names) const
		{
		  return names.end() != std::find(names.begin(), names.end(), u->Name());
		}
	};

  // predicate to determine if the given atom or molecule has the given name
  template <class U>
	class IsName : public std::binary_function<U, std::string, bool> {
	  public:
		bool operator() (const U u, const std::string name) const {
		  return u->Name() == name; 
		}
	};

  // returns an iterator to the first occurence of a member with the given name
  template <typename Iter>
	Iter FindByName (Iter first, Iter last, const std::string& name) {
	  typedef typename std::iterator_traits<Iter>::value_type val_t;
	  return std::find_if (first, last, std::bind2nd(IsName<val_t>(), name));
	}

  template <class U> 
	void KeepByName (U& u, std::string& name) {
	  u.erase(
		  remove_if(u.begin(), u.end(), not1(std::bind2nd(IsName<typename U::value_type>(), name))), u.end()
		  );
	  return;
	}

  // keep elements of the vector with names matching one of those in the list
  template <class U> 
	void KeepByNames (U& u, std::vector<std::string>& names) {
	  u.erase(
		  remove_if(u.begin(), u.end(), not1(std::bind2nd(NameInList<typename U::value_type>(), names))), u.end());
	  return;
	}

  // remove all elements that have the given name
  template <class U> 
	void RemoveByName (U& u, std::string& name) {
	  u.erase(
		  remove_if(u.begin(), u.end(), std::bind2nd(IsName<typename U::value_type>(), name)), u.end()
		  );
	  return;
	}

  // remove all elements that have names matching any of those in the list of names supplied
  template <class U> 
	void RemoveByNames (U& u, std::vector<std::string>& names) {
	  u.erase(
		  remove_if(u.begin(), u.end(), std::bind2nd(NameInList<typename U::value_type>(), names)), u.end());
	  return;
	}

}	// namespace md_utility


namespace histogram {

  template <typename T>
	class ValueBetween : public std::binary_function<T,std::pair<T,T>,bool> {
	  public:
		bool operator() (const T& val, const std::pair<T,T>& test) const {
		  return (val > test.first) && (val < test.second);
		}
	}; // value between - predicate

  template <typename Iter> 
	//std::vector< std::pair<typename std::iterator_traits<Iter>::value_type, int> > 
	std::vector< std::pair<double, int> > 
	Histogram (Iter first, Iter last, const int num_bins) {

	  typedef typename std::iterator_traits<Iter>::value_type val_t;

	  // sorts the data - note this is a destructive operation!
	  std::sort (first, last);

	  // instead of requiring a min/max value to be supplied, we just take the smallest and highest values in the data set
	  val_t max = *std::max_element(first,last);
	  val_t min = *std::min_element(first,last);
	  val_t bin_size = (max - min)/((val_t)num_bins);

	  std::vector< std::pair<val_t,int> > histogram;

	  std::pair<val_t,val_t> test;
	  std::pair<val_t,int> result;

	  for (val_t bin = min; bin < max; bin += bin_size) {
		test.first = bin;
		test.second = bin + bin_size;

		result.first = bin;
		result.second = (int) std::count_if (first, last, std::bind2nd(ValueBetween<double>(), test));

		histogram.push_back(result);
	  }

	  return histogram;
	}	// Histogram converter





  /*******************************************************************************************************************************************************/
  /***************************************** Histograms **************************************************************************************************/


  /* 1-d histogram functor */
  template <class T>
	class Histogram1D : public std::unary_function<T,bool>
  {
	private:
	  typedef double histo_element_t;

	  T _min, _max, _res;
	  int _size;
	  double _access_count;
	  std::vector<histo_element_t> _histogram;

	  // Returns the bin for a given value
	  int Bin (const T t) const { return int ((t - _min)/_res); }
	  bool InBounds (const T t) const { return (t >= _min && t <= _max); }

	public:

	  Histogram1D (const T min, const T max, const T res) 
		: _min(min), _max(max), _res(res), 
		_size(int((max - min)/res) + 1), _access_count(0.0)
	{ 
	  if (min > max) {
		printf ("Check the limits given to Histogram1D - minimum is greater than the maximum!!\n");
		exit(1);
	  }
	  _histogram.resize(_size, histo_element_t(0));
	}

	  bool operator() (const T t) {
		bool ret = false;
		if (InBounds(t)) 
		{
		  _histogram [Bin(t)]++;
		  _access_count++;
		  if (_access_count < 0.0) 
		  {
			printf ("weird value started with: %f in the histogram operator()\n", t);
			std::cout << "bin was " << Bin(t) << " with a population of " << _histogram[Bin(t)] << std::endl;
			std::cout << "the _access_count is: " << _access_count << std::endl;
			exit(1);
		  }
		  ret = true;
		}
		return ret;
	  }

	  double Count () const { return _access_count; }
	  int Size () const { return _size; }
	  T Max () const { return _max; }
	  T Min () const { return _min; }
	  T Resolution () const { return _res; }

	  // returns the population of a single bin given a value
	  histo_element_t Population (const T t) const { 
		if (!InBounds(t)) {
		  printf ("Population requested for a value outside of the histogram limits\n");
		  exit(1);
		}
		return _histogram[this->Bin(t)];
	  }
  };	// 1D Histogram


  /* A 2-dimensional histogram functor. Each application of the histogram will bin the given values. Initial setup requires some parameters for the histogram (size, resolution, etc) */
  template <class T>
	class Histogram2D : public std::binary_function<T,T,bool>
  {
	public:

	  typedef std::pair<T,T>	pair_t;
	  pair_t min;						// minimum values for each dimension
	  pair_t max;						// maximum value the histogram can bin in each dimension
	  pair_t resolution;				// resolution for each dimension
	  std::pair<int,int> size;					// dimensions of the 2-d data (number of histogram bins)
	  std::vector<double> counts;					// A running total of each time a bin was updated in the 1st dimension of the histogram i.e. access count

	  // Initialization with [min, max, resolution]
	  Histogram2D (const pair_t& minima, const pair_t& maxima, const pair_t& resolutions) 
		: min(minima), max(maxima), resolution(resolutions)
	  {
		size.first = int((max.first - min.first)/resolution.first) + 1;
		size.second = int((max.second - min.second)/resolution.second) + 1;

		/* 
		   printf ("size = <%d,%d>\nmax = <%f,%f\n,min = <%f,%f>\nres = <%f,%f>\n",
		   size.first, size.second, max.first, max.second, min.first, min.second, resolution.first, resolution.second);
		 */

		_histogram.clear();
		_histogram.resize (size.first, Histogram_t (size.second, 0));

		counts.resize(size.first, 0);
	  }

	  void operator() (const T& a, const T& b) {
		if (CheckLimits(a,b))
		{
		  // find which bin is to be updated
		  bins new_bins = Bin(a,b);
		  // Update the correct bin
		  _histogram[new_bins.first][new_bins.second]++;
		  // update the number of times that this particular 1st-dimensions has been accessed (total # of bin updates)
		  counts[new_bins.first]++;
		}
		return;
	  }

	  /* Increment the bins given by the two indices without doing any checks */
	  void Shove (const int a, const int b)
	  {
		_histogram[a][b]++;
		counts[a]++;
		return;
	  }
	  // Return the element of the histogram
	  int Element (const int x, const int y) const { return _histogram[x][y]; }

	  // returns the population of a single bin given a value
	  double Population (const T& a, const T& b) const 
	  { 
		bins temp = Bin(a, b);
		return _histogram[temp.first][temp.second]; 
	  }		

	  //int Count (const int i) const { return counts[i]; }
	  double Count (const T& i) const { return counts[(i-min.first)/resolution.first]; }

	private:
	  typedef std::vector<double> Histogram_t;
	  std::vector<Histogram_t> _histogram;		// 2-d container/histogram

	  typedef unsigned int bin;
	  typedef std::pair<bin,bin> bins;

	  /* calculates the bins into which a value will be placed */
	  bins Bin (const T& a, const T& b) const
	  { 
		return std::make_pair( bin((a-min.first)/resolution.first), bin((b-min.second)/resolution.second) );
	  }

	  /* checks if the given value pair is within the limits of the histogram */
	  bool CheckLimits (const T& a, const T& b) const
	  { 
		return a > min.first && a < max.first && b > min.second && b < max.second;
	  }

  };	// histogram 2D

}	// namespace histogram

#endif
