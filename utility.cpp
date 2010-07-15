#include "utility.h"

namespace md_utility {


  /*******************************************************************************************************************************************************/
  /***************************************** Dealing with Pairs ******************************************************************************************/


  template <class Iter>
	Iter PairListMember (const typename std::iterator_traits<Iter>::value_type& p, Iter first, Iter last)
	{
	  return find_if (first, last, std::bind1st(EqualPairs<typename std::iterator_traits<Iter>::value_type>(), p));
	}

  template <class U>
	void KeepByName (U& u, std::string& name) {
	  u.erase(
		  remove_if(u.begin(), u.end(), not1(std::bind2nd(IsName<typename U::value_type>(), name))), u.end()
		  );
	  return;
	}

  template <class U>
	void KeepByNames (U& u, std::vector<std::string>& names) {
	  u.erase(
		  remove_if(u.begin(), u.end(), not1(std::bind2nd(md_utilities::NameInList<typename U::value_type>(), names))), u.end());
	  return;
	}


  template <class U>
	void RemoveByName (U& u, std::string& name) {
	  u.erase(
		  remove_if(u.begin(), u.end(), std::bind2nd(IsName<typename U::value_type>(), name)), u.end()
		  );
	  return;
	}


  template <class U>
	void RemoveByNames (U& u, std::vector<std::string>& names) {
	  u.erase(
		  remove_if(u.begin(), u.end(), std::bind2nd(md_utilities::NameInList<typename U::value_type>(), names)), u.end());
	  return;
	}

} // namespace md_utilities



template <typename Iter> std::vector< std::pair<double,int> > 
//std::vector< std::pair<typename std::iterator_traits<Iter>::value_type, int> > 
histogram::Histogram (Iter first, Iter last, const int num_bins) {

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






