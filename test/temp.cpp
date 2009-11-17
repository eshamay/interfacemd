#include <iostream>
#include "../utility.h"

using namespace std;

class A {
  public:
  A (string n) : name(n), id(++A::numA) { return; }
  ~A () { return; }

  static int numA;
  string name;
  int id;
};

int A::numA = 0;

MAKE_PREDICATE(test_h2o, A *, t->id >= 3);
MAKE_FUNCTOR(change_name, void, A *, if (t->name == "h2o") t->name = "blat";)
MAKE_FUNCTOR(print_name, void, A *, printf ("%d %s\n", t->id, t->name.c_str());)
double add_three (double i) { return i+3.0; }

int main (int argc, char **argv) {

  std::vector<A *> a;
  a.push_back(new A("h2o"));
  a.push_back(new A("foo"));
  a.push_back(new A("bar"));
  a.push_back(new A("foo"));
  a.push_back(new A("h2o"));
  a.push_back(new A("baz"));
  a.push_back(new A("h2o"));
  a.push_back(new A("bar"));

FOR_EACH(a,print_name());
FOR_EACH(a,change_name());
FOR_EACH(a,print_name());

return 0;
}

