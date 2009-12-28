#include <iomanip>
#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

typedef std::ostream_iterator<double> oi_t;
typedef std::istream_iterator<double> ii_t;

struct RandomGenerator
{
  double operator() ()
  {
	return (double)rand()/RAND_MAX;
  }
};

int main () {

  // generate some data to test
  std::vector<double> vd(20);
  generate(vd.begin(), vd.end(), RandomGenerator());


  // perform output to a binary file
  ofstream output ("temp.bin", ios::binary);
  oi_t oi (output, " ");
  output << setprecision(9);
  copy (vd.begin(), vd.end(), oi_t(std::cout, "\n"));

  output.close();

  // input from the binary file to a container
  ifstream input ("temp.bin", ios::binary);
  ii_t ii(input);
  ii_t ii_end;
  std::vector<double> data(ii, ii_end);
  input.close();

  // output data to screen to verify/compare the results
  for (int i = 0; i < vd.size(); i++)
    printf ("%8.4f  %8.4f\n", vd[i], data[i]);

  printf ("vd.size() = %d\tvi.size() = %d\n", vd.size(), data.size());
  return 0;
}

