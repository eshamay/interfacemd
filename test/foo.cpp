#include "../utility.h"

using namespace std;

int main () {

  typedef std::pair<double,int> pair_t;
  typedef std::vector<pair_t> vp_t;

  vp_t vp;

  for (int i = 0; i < 10; i++) {
	vp.push_back(std::make_pair((double)i*2.4, i));
  }

  std::vector<double> vi(10);

  std::transform (vp.begin(), vp.end(), vi.begin(), pair_utility::pair_ref_1st<pair_t>());

  for (vp_t::iterator it = vp.begin(); it != vp.end(); it++) {
	printf ("%f %d\n", it->first, it->second);
  }

  std::copy (vi.begin(), vi.end(), std::ostream_iterator<double>(std::cout, " "));

  return 0;
}
