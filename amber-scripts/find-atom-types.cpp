#include "../pdbsystem.h"

using namespace std;

int main (int argc, char **argv) {

	if (argc < 3) {
		printf ("Command : <program> <pdb-file-path> <output-choice>\n");
		printf ("Where output-choice is either of the following:\n");
		printf ("\t1) Atom table\n");
		printf ("\t2) Perturbed atom table\n");
		exit(1);
	}

	md_files::PDBSystem pdb (argv[1]);

	int output_choice = atoi(argv[2]);

	if (output_choice == 1)
		printf ("!entry.qtz.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg\n");
	else
		printf ("!entry.h2o.unit.atomspertinfo table  str pname  str ptype  int ptypex  int pelmnt  dbl pchg\n");

	int i = 0;
	for (Atom_it it = pdb.begin(); it != pdb.end(); it++) {

		(*it)->SetAtomProperties();

		VecR pos = (*it)->Position();
		std::string atom_type;
		int num;
		double charge;

		if ((*it)->Element() == Atom::H) {
			atom_type = "HC";
			charge =  0.06;
			num = 1;
		}

		if ((*it)->Element() == Atom::C) {
			atom_type = "CT";
			charge =  -0.12;
			num = 12;
		}

		if ((*it)->Element() == Atom::O)
		{
			// check for OS or OH
			if (pos[y] > -8.5)
			{
				atom_type = "OH";
				charge = -0.8715;
			}
			else
			{
				atom_type = "OS";
				charge = -0.8715;
			}
			num = 8;
		}

		if ((*it)->Element() == Atom::Si)
		{
			atom_type = "SI";
			charge = 1.428;
			num = 14;
		}

		if (output_choice == 1)
			printf (" \"%s\" \"%s\" 0 1 131075 %d %d %8.6f\n", (*it)->Name().c_str(), atom_type.c_str(), i+1, num, charge);
		if (output_choice == 2)
			printf (" \"%s\" \"%s\" 0 -1 0.0\n", (*it)->Name().c_str(), atom_type.c_str());

		i++;
	}

	return 0;
}
