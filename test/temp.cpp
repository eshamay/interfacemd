#include <iostream>
#include "../utility.h"

enum names {"a", "b", "c", "d", "e"};

void print_name (int a) {
  printf ("%s\n", names(a));
  return;
}

int main () {

  vector<int> a;
  for (int i = 0; i < 20; i++) {
	a.push_back (rand() % 5);
  }
  FOR_EACH(a, print_name);

  return 0;
}
