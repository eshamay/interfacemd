#include "../pdbfile.h"

using namespace std;

int main () {

  PDBFile pdb ("quartz_100-slab+H.pdb");

  Molecule * mol = pdb.Molecules(0);
  std::vector<Atom *> atoms = mol->Atoms();

  printf ("!entry.qtz.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg\n");
  for (int i = 0; i < atoms.size(); i++)
  {
	Atom * atom = atoms[i];
	VecR pos = atom->Position();

	std::string atom_type, atom_name;
	atom_name = atom->Name();
	double charge;
	int num;

	const char * a = atom_name.c_str();
	if (a[0] == 'H')
	{
	  if (pos[y] > -8.5)
	  { 
		atom_type = "HO"; 
		charge = 0.41;
		//charge = 0.32;
	  }
	  else 
	  {
		atom_type = "HA";
		charge = 0.41;
		//charge = -0.15;
	  }
	  num = 1;
	}

	if (a[0] == 'O')
	{
	  if (pos[y] > -8.5)
	  {
		atom_type = "OH";
		charge = -1.07;
		//charge = -.54;
	  }
	  else
	  {
		atom_type = "OS";
		charge = -.66;
		//charge = -.5425;
	  }
	  num = 8;
	}

	if (a[0] == 'S')
	{
		if (pos[y] < -16.0)
		{
		  atom_type = "SH";
		  //charge = 0.8;
		  charge = 0.5;
		}
		else {
		  atom_type = "SI";
		  //charge = 1.08;
		  charge = 1.32;
		}
		num = 14;
	}
		  
	printf (" \"%s\" \"%s\" 0 1 131075 %d %d %8.6f\n", atom_name.c_str(), atom_type.c_str(), i+1, num, charge);
	//printf (" \"%s\" \"%s\" 0 -1 0.0\n", atom_name.c_str(), atom_type.c_str());
  }

  return 0;
}
