#include "test.h"

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


Tester::Tester (WaterSystemParams& wsp)
  :	
	Analyzer<XYZSystem> (wsp)
{
  return;
}

void Tester::Setup () {

  return;
}

void Tester::Analysis () {

  LoadAll();
  //std::string name ("so2");
  //KeepByName<Mol_ptr_vec> (int_mols, name);

  for (Mol_it it = sys_mols.begin(); it != sys_mols.end(); it++) {
	sys->CalcDipole(*it);
	dipoles.push_back((*it)->Dipole().Magnitude());
  }

  return;

}

void Tester::DataOutput (const unsigned int timestep) { 

  typedef std::vector< std::pair<double, int> > histo_t;
  typedef histo_t::const_iterator histo_it;

  histo_t histo = histogram::Histogram (dipoles.begin(), dipoles.end(), 100);

  for (histo_it it = histo.begin(); it != histo.end(); it++) {
	printf ("% 13.3f % 13d\n", it->first, it->second);
  }

  return; 
}



int main () {

  libconfig::Config cfg;
  cfg.readFile("system.cfg");

  std::string filename = cfg.lookup("analysis.test.filename");
  libconfig::Setting &analysis = cfg.lookup("analysis");
  analysis.add("filename", libconfig::Setting::TypeString) = filename;

  WaterSystemParams wsp (cfg);

  Tester test (wsp);

  test.SystemAnalysis ();

  return 0;
}

