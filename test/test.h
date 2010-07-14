#ifndef TEST_H_
#define TEST_H_

#include "../analysis.h"
#include "../utility.h"

class Tester : public Analyzer<XYZSystem>
{

  public:
	Tester (WaterSystemParams& params);

  private:

	void Setup ();
	void Analysis ();
	void DataOutput (const unsigned int timestep);
	void PostAnalysis () { return; }

	std::vector<double>	dipoles;

};

namespace histogram {

  // predicate test if a value falls between a low/high pair
  template <typename T>
	class ValueBetween : public std::binary_function<T,std::pair<T,T>,bool> {
	  public:
		bool operator() (const T& val, const std::pair<T,T>& test) const {
		  return (val > test.first) && (val < test.second);
		}
	}; // value between - predicate


  // Returns a histogram of the given scalar array
  // return[i].first == bin
  // return[i].second == population
  template <typename Iter> 
	//std::vector< std::pair<typename std::iterator_traits<Iter>::value_type, int> > 
	std::vector< std::pair<double, int> > 
	Histogram (Iter first, Iter last, const int num_bins);

}	// namespace histogram

#endif

