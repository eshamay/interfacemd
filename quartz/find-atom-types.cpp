#include "../pdbfile.h"

using namespace std;

int main () {

	PDBFile pdb ("octadecane.pdb");

	Molecule * mol = pdb.Molecules(0);
	Atom_ptr_vec atoms = mol->Atoms();

	printf ("!entry.qtz.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg\n");
	int i = 0;
	for (Atom_it atom_i = mol->begin(); atom_i != mol->end(); atom_i++)
	{
		Atom * atom = (*atom_i);
		VecR pos = atom->Position();

		std::string atom_type, atom_name;
		atom_name = atom->Name();
		double charge;
		int num;

		const char * a = atom_name.c_str();
		if (a[0] == 'H') {
			atom_type = "HC";
			charge =  0.06;
			num = 1;
		}

		if (a[0] == 'C') {
			atom_type = "CT";
			charge =  -0.12;
			num = 12;
		}

		if (a[0] == 'O')
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

		if (a[0] == 'S')
		{
			atom_type = "SI";
			charge = 1.428;
			num = 14;
		}

		//printf (" \"%s\" \"%s\" 0 1 131075 %d %d %8.6f\n", atom_name.c_str(), atom_type.c_str(), i+1, num, charge);
		printf (" \"%s\" \"%s\" 0 -1 0.0\n", atom_name.c_str(), atom_type.c_str());
		i++;
	}

	return 0;
}
