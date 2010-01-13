#ifndef __UTIL_H
#define __UTIL_H

#include <vector>
#include <algorithm>
#include <utility>
#include <functional>
#include <iterator>
#include <iostream>

/*
void DEBUGMSG (const std::string msg)
{
#ifdef __DEBUG__
  std::cout << msg << std::endl;
#endif
  return;
}
*/

template <class Iter>
Iter PairListMember (const typename std::iterator_traits<Iter>::value_type& p, Iter first, Iter last);
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

template<class C>
struct Inserter {
  std::back_insert_iterator<C> in;
  Inserter(C& c) : in(c) {}
  void operator()(const std::pair<typename C::value_type, typename C::value_type>& p)
  {
	*in++ = p.first;
	*in++ = p.second;
  }
};

  template<class C>
Inserter<C> make_inserter(C& c)
{ 
  return Inserter<C>(c); 
}

/* Takes a container of pairs (pairlist) and flattens it such that each pair element is now an element of the new list */
#define FLATTEN_PAIRLIST(pairlist,result)	\
  std::for_each(pairlist.begin(), pairlist.end(), make_inserter(result));

/* sorts and finds unique elements of a given list and stores them into the result */
#define UNIQUE_PAIRLIST_ELEMENTS(pairlist,result)	\
  FLATTEN_PAIRLIST (pairlist, result)	\
  sort(result.begin(), result.end());	\
  result.resize( unique(result.begin(), result.end()) - result.begin());




/*******************************************************************************************************************************************************/
/***************************************** Histograms **************************************************************************************************/

/* 1-d histogram functor */
template <class T>
class Histogram1D : public std::unary_function<T,bool>
{
  private:
	typedef T histo_element_t;

	T _min, _max, _res;
	int _size;
	int _access_count;
	std::vector<histo_element_t> _histogram;

	int Bin (const T t) const { return int ((t - _min)/_res); }

  public:

	Histogram1D (const T min, const T max, const T res) : _min(min), _max(max), _res(res), 
														  _size(int((max - min)/res) + 1),
														  _access_count(0)
  { 
	if (min > max) {
	  printf ("Check the limits given to Histogram1D - minimum is greater than the maximum!!\n");
	  exit(1);
	}
	_histogram.resize(_size, histo_element_t(0));
  }


	bool operator() (const T t) {
	  bool ret = false;
	  if (t < _min || t > _max) { 
		ret = false;				// out of bounds issue...
	  }
	  else {
		_histogram [Bin(t)]++;
		_access_count++;
		ret = true;
	  }

	  return ret;
	}

	int Count () const { return _access_count; }
	int Size () const { return _size; }
	int Max () const { return _max; }
	int Min () const { return _min; }
	int Resolution () const { return _res; }

	int Population (const T t) const { return _histogram[this->Bin(t)]; }		// returns the population of a single bin given a value
	std::vector<histo_element_t>& Histogram () { return _histogram; }

};
	

/* A 2-dimensional histogram functor. Each application of the histogram will bin the given values. Initial setup requires some parameters for the histogram (size, resolution, etc) */
template <class T>
class Histogram2D : public std::binary_function<T,T,bool>
{
  public:

  typedef std::pair<T,T>	pair_t;
  pair_t max;						// maximum value the histogram can bin in each dimension
  pair_t min;						// minimum values for each dimension
  pair_t resolution;				// resolution for each dimension
  std::pair<int,int> size;					// dimensions of the 2-d data (number of histogram bins)
  std::vector<int> counts;					// A running total of each time a bin was updated in the 1st dimension of the histogram i.e. access count

  // Initialization
  Histogram2D (const pair_t& maxima, const pair_t& resolutions, const pair_t& minima = std::make_pair(T(0), T(0))) 
	: max(maxima), min(minima), resolution(resolutions)
  {
	size.first = int((max.first - min.first)/resolution.first) + 1;
	size.second = int((max.second - min.second)/resolution.second) + 1;

	/* 
	printf ("size = <%d,%d>\nmax = <%f,%f\n,min = <%f,%f>\nres = <%f,%f>\n",
		size.first, size.second, max.first, max.second, min.first, min.second, resolution.first, resolution.second);
	*/

	histogram.clear();
	histogram.resize (size.first, Histogram_t (size.second, 0));

	counts.resize(size.first, 0);
  }

  void operator() (const T& a, const T& b) {
	if (CheckLimits(a,b))
	{
	  // find which bin is to be updated
	  bins new_bins = Bin(a,b);
	  // Update the correct bin
	  histogram[new_bins.first][new_bins.second]++;
	  // update the number of times that this particular 1st-dimensions has been accessed (total # of bin updates)
	  counts[new_bins.first]++;
	}
	return;
  }

  /* Increment the bins given by the two indices without doing any checks */
  void Shove (const int a, const int b)
  {
	histogram[a][b]++;
	counts[a]++;
	return;
  }
  // Return the element of the histogram
  int Element (const int x, const int y) const { return histogram[x][y]; }
  int Count (const int i) const { return counts[i]; }

  private:
  typedef std::vector<int> Histogram_t;
  std::vector<Histogram_t> histogram;		// 2-d container/histogram
 
  typedef unsigned int bin;
  typedef std::pair<bin,bin> bins;

  /* calculates the bins into which a value will be placed */
  bins Bin (const T& a, const T& b)
  { 
	return std::make_pair( bin((a-min.first)/resolution.first), bin((b-min.second)/resolution.second) );
  }

  /* checks if the given value pair is within the limits of the histogram */
  bool CheckLimits (const T& a, const T& b) const
  { 
	return a > min.first && a < max.first && b > min.second && b < max.second;
  }

};

/*******************************************************************************************************************************************************/
/***************************************** Simplifiers *************************************************************************************************/

// These RUN macros should clean up code for checking through vectors and arrays
#define RUN(list)   \
  for (size_t i = 0; i < list.size(); i++)

#define RUN2(list)  \
  for (size_t j = 0; j < list.size(); j++)

// A destructive remove_if that erases the given list with a new list containing only the elements that don't match the predicate op
#define D_REMOVE_IF(vec,pred)	\
  vec.erase (std::remove_if (vec.begin(), vec.end(), pred), vec.end() );

// A general method for making new functors
#define MAKE_FUNCTOR(name,return_type,arg_type,code)	\
  class name {	\
	public:	\
	  return_type operator() (arg_type t) { code }	\
  };

// Creates a new predicate object that has a name as given, and takes a value of the given type (can be an object), and returns the value of the 'code' provided. code must return a boolean for this to work.
#define MAKE_PREDICATE(name,arg_type,code)	\
  MAKE_FUNCTOR(name, bool, arg_type, return code;);

#define FOR_EACH(list,fn)	\
  std::for_each (list.begin(), list.end(), fn);

// Maps OP onto A and stores the result into B
#define MAP_TO(A,B,OP)	\
  std::transform(A.begin(), A.end(), B.begin(), OP);

// Maps OP onto A and stores the result back into A - destructive function
#define MAP(A,OP)	\
  MAP_TO(A,A,OP);


/*******************************************************************************************************************************************************/
/***************************************** Binary file I/O functors/iterators *************************************************************************************************/

/*
// output-stream iterator for writing binary data to an output
template<typename T>
struct oi_t: public std::iterator<output_iterator_tag, void, void, void, void>
{
  oi_t(std::ostream& str)
    :m_str(str)
  {}
  oi_t operator++()   {return *this;}  // increment does not do anything.
  oi_t operator++(int){return *this;}
  oi_t operator*()    {return *this;}  // Dereference returns a reference to this
                                       // So that when the assignment is done we
                                       // actually write the data from this class
  oi_t operator=(T const& data)
  {
    // Write the data in a binary format
    m_str.write(reinterpret_cast<char const*>(&data),sizeof(T));
  }

  private:
    std::ostream&   m_str;
};

// istream iterator for reading in data from a binary file
template<typename T>
struct ii_t: public iterator<input_iterator_tag, void, void, void, void>
{
  ii_t(std::istream& str)
    :m_str(&str)
  {}
  ii_t()
    :m_str(NULL)
  {}
  ii_t operator++()   {return *this;}  // increment does nothing.
  ii_t operator++(int){return *this;}
  T& operator*()
  {
    // On the de-reference we actuall read the data into a local //// static ////
    // Thus we can return a reference
    static T result;
    m_str->read(reinterpret_cast<char*>(&result),sizeof(T));
    return result;
  }
  // If either iterator has a NULL pointer then it is the end() of stream iterator.
  // Input iterators are only equal if they have read past the end of stream.
  bool operator!=(ii_t const& rhs)
  {
      bool lhsPastEnd = (m_str == NULL)     || (!m_str->good());
      bool rhsPastEnd = (rhs.m_str == NULL) || (!rhs.m_str->good());

      return !(lhsPastEnd && rhsPastEnd);
  } 

  private:
    std::istream*   m_str;
};
*/

#endif
