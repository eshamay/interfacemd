#include "mpitest.h"

int main (int argc, char **argv) {

	MPIMolSystem sys (&argc, &argv, PRMTOP, MDCRD, FORCE);

		std::cout << sys.ID() << std::endl;
		sys.Molecules(0)->Print();
		std::cout << sys.ID() << std::endl;
		sys.Molecules(800)->Print();

return 0;
}
