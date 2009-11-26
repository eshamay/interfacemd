#include <iostream>
#include "../analysis.h"

using namespace std;

typedef std::pair<int,int> int_pair;
typedef std::map<int_pair, int>	int_map;

int main () {

  //WaterSystemParams wsp ("temp.out");
  //Analyzer an(wsp);


  int_map a;

  vector<int_pair> b;
  b.push_back (std::make_pair(1,2));
  b.push_back (std::make_pair(2,2));
  b.push_back (std::make_pair(3,3));
  b.push_back (std::make_pair(4,2));
  b.push_back (std::make_pair(5,3));

  for (vector<int_pair>::iterator it = b.begin(); it != b.end(); it++)
  {
	a[*it] = it->first + it->second;
  }

  a[std::make_pair(4,2)] = 7;

  for (int_map::iterator it = a.begin(); it != a.end(); it++) {
	std::cout << it->first.first << " + " << it->first.second << " = " << it->second << std::endl;
  }

  return 0;
}
