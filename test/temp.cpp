#include <iostream>
#include <vector>
#include <functional>
#include <utility>
#include "../utility.h"

using namespace std;

typedef pair<int,int> pair_t;

int main () {

  vector< pair_t > a;
  a.push_back (make_pair(1,2));
  a.push_back (make_pair(2,2));
  a.push_back (make_pair(3,2));
  a.push_back (make_pair(1,3));
  a.push_back (make_pair(3,3));
  a.push_back (make_pair(4,2));

  vector< pair_t > b;
  b.push_back (make_pair(5,1));

  bool in_list = PAIR_IN_LIST(make_pair(3,1), a);

  if (in_list) 
	printf ("found it\n");
  else 
	printf ("didn't find it\n");

  return 0;
}
