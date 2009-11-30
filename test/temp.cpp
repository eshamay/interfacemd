#include <iostream>
#include "../analysis.h"

using namespace std;


/* a predicate for name sorting */
template <class T>
struct Name_equal_pred {
    bool operator()(const T &left, const T &right) {
        return left->Name() < right->Name();
    }
};


int main () {

  WaterSystemParams wsp ("temp.out");
  Analyzer an(wsp);

  an.LoadAll();
  Atom_ptr_vec atoms = an.Atoms();

  std::vector<string> names;
  names.push_back("O");
  names.push_back("C");
  names.push_back("H1");
  names.push_back("H2");
  names.push_back("Cl1");
  names.push_back("Cl3");
  names.push_back("Cl4");

  KEEP_BY_NAMES(atoms, Atom *, names);

  RUN (atoms)
	std::cout << atoms[i]->Name() << std::endl;

  return 0;
}
