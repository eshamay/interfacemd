#include "gmxtest.h"

int main () {

	GMXSystem sys ("gro", "trr");

	sys.Molecules(0)->Print();
	sys.LoadNext();
	sys.Molecules(0)->Print();

  return 0;

}
