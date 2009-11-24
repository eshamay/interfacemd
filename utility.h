#ifndef __UTIL_H_
#define __UTIL_H_

#include <vector>
#include <algorithm>
#include <utility>
#include <functional>


/*******************************************************************************************************************************************************/
/***************************************** Dealing with Pairs ******************************************************************************************/


/* A functor that takes a std::pair as a constructor argument, and all applications of the tester will test pairs against the initial one given. */
template <class T>
struct EqualPairs : public std::binary_function<T,T,bool> 
{
  bool operator() (const T& lhs, const T& rhs) const 
  {
	return lhs == rhs || lhs.first == rhs.second && lhs.second == rhs.first;
  }
};

/* Searches for the supplied pair p in the sequence of pairs from first to last */
template <class Iter>
Iter PairListMember (const typename std::iterator_traits<Iter>::value_type& p, Iter first, Iter last)
{
  return find_if (first, last, std::bind1st(EqualPairs<typename std::iterator_traits<Iter>::value_type>(), p));
}

/* A macro to simplify the PairListMember function to act as a predicate for member existence in a list */
#define PAIR_IN_LIST(pair,list)	\
  list.end() != PairListMember(pair, list.begin(), list.end())

#define PAIR_LIST_MEMBER



/*******************************************************************************************************************************************************/
/***************************************** Histograms **************************************************************************************************/

/* A 2-dimensional histogram functor. Each application of the histogram will bin the given values. Initial setup requires some parameters for the histogram (size, resolution, etc) */
template <class T>
struct 2DHistogram : public std::binary_function<T,T,bool>
{
  std::pair<T,T> max;						// maximum value the histogram can bin in each dimension
  std::pair<T,T> min;						// minimum values for each dimension
  std::pair<T,T> resolution;				// resolution for each dimension
  std::pair<int,int> size;					// dimensions of the 2-d data (number of histogram bins)

  typedef std::vector<T> Histogram_t;
  std::vector<Histogram> histogram;			// 2-d container/histogram

  typedef int bin;
  typedef std::pair<bin,bin> bins;

  // Initialization
  2DHistogram (const std::pair<T,T>& maxima, const std::pair<T,T>& resolutions, const std::pair<T,T>& minima = make_pair(T(0), T(0))) 
	: max(maxima), resolution(resolutions), min(minima)
  {
	size.first = int((max.first - min.first)/resolution.first) + 1;
	size.second = int((max.second - min.second)/resolution.second) + 1;

	histogram.clear();
	histogram.resize (size.first, Histogram_t (size.second, 0));
  }

  /* calculates the bins into which a value will be placed */
  bins Bin (const T& a, const T& b)
  { return make_pair( bin((a-min.first)/resolution.first), bin((b-min.second)/resolution.second) ); }

  /* checks if the given value pair is within the limits of the histogram */
  bool CheckLimits (const T& a, const T& b) const
  { return a > min.first && a < max.first && b > min.second && b < max.second; }

  bool operator() (const T& a, const T& b) {
	if (CheckLimits(a,b))
	{
	  bins b = Bin(a,b);
	  histogram[b.first][b.second]++;
	  return true;
	}
	else { return false; }
  }
};

/*******************************************************************************************************************************************************/
/***************************************** Simplifiers *************************************************************************************************/

// These RUN macros should clean up code for checking through vectors and arrays
#define RUN(list)   \
  for (size_t i = 0; i < list.size(); i++)

#define RUN2(list)  \
  for (size_t j = 0; j < list.size(); j++)


#define REMOVE_IF(vec,pred)	\
  vec.erase (std::remove_if (vec.begin(), vec.end(), pred()), vec.end() );

// A general method for making new functors
#define MAKE_FUNCTOR(name,return_type,arg_type,code)	\
  class name {	\
	public:	\
	  return_type operator() (arg_type t) { code }	\
  };

// Creates a new predicate object that has a name as given, and takes a value of the given type (can be an object), and returns the value of the 'code' provided. code must return a boolean for this to work.
#define MAKE_PREDICATE(name,arg_type,code)	\
  MAKE_FUNCTOR(name, bool, arg_type, return code;);

#define FOR_EACH(name,fn)	\
  std::for_each (name.begin(), name.end(), fn);

// Maps OP onto A and stores the result into B
#define MAP_TO(A,B,OP)	\
  std::transform(A.begin(), A.end(), B.begin(), OP);

// Maps OP onto A and stores the result back into A - destructive function
#define MAP(A,OP)	\
  MAP_TO(A,A,OP);

#endif
