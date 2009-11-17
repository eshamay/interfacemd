#include <iostream>
#include "../utility.h"

using namespace std;

class Thing {
  public:
  Thing (string n) : name(n), id(++Thing::numThing) { return; }
  ~Thing () { return; }

  static int numThing;
  string name;
  int id;
};

int Thing::numThing = 0;

// A way to change the name of an object
MAKE_FUNCTOR(change_name, void, Thing *, if (t->name == "h2o") t->name = "blat";)

// A quick way to print all the names in a vector
MAKE_FUNCTOR(print_thing, void, Thing *, printf ("%d %s\n", t->id, t->name.c_str());)

// A predicate functor for testing if the argument is a "water"
MAKE_PREDICATE(bar_p, Thing *, t->name == "bar")

// A functor to change the id of all the 'things'
MAKE_FUNCTOR(change_id, void, Thing *, t->id += 10;)

int main (int argc, char **argv) {

  // create a vector of type Thing * to hold a bunch of 'things'
  std::vector<Thing *> a;
  a.push_back(new Thing("h2o"));
  a.push_back(new Thing("foo"));
  a.push_back(new Thing("bar"));
  a.push_back(new Thing("foo"));
  a.push_back(new Thing("h2o"));
  a.push_back(new Thing("baz"));
  a.push_back(new Thing("h2o"));
  a.push_back(new Thing("bar"));

 
  // Print out all the names of the things
  FOR_EACH(a,print_thing())
  // Change the names of "h2o" things to "blat"
  FOR_EACH(a,change_name())
// Print out the new names
  FOR_EACH(a,print_thing())
  // Remove all the "bar" objects
  REMOVE_IF(a, bar_p)
  // Print out the remaining names
  FOR_EACH(a, print_thing())
  // Change all the id's by adding 10
  FOR_EACH(a,change_id())
  // Print out the new ids
  FOR_EACH(a, print_thing())

return 0;
}

