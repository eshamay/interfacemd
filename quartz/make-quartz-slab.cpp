#include "../pdbfile.h"
#include "quartz-utils.h"


using namespace std;

int main () {
  PDBFile pdb ("quartz_100_short.pdb");

  // create the quartz unitcell
  Molecule * uc = pdb.Molecules(0);

  //make copies of the quartz unitcell and shift them around
  VecR x_shift (4.91, 0.0, 0.0);
  VecR y_shift (0.0, 0.0, 0.0);
  VecR z_shift (0.0, 5.402, 0.0);

  std::vector<Molecule *> mols;

  VecR yshift, xshift, zshift;
  for (int i = 0; i < 1; i++) {
	xshift.Zero();
	zshift.Zero();

	for (int j = 0; j < 6; j++) {
	  zshift.Zero();

	  for (int k = 0; k < 6; k++) {

		// make a copy of the unitcell
		Molecule * copy = new Molecule();
		copy->Name("qtz");
		// make copies of all the atoms in the unitcell
		for (int atom = 0; atom < uc->size(); atom++) {
		  Atom * pa = new Atom (*uc->Atoms(atom));
		  copy->AddAtom(pa);
		}

		// shift the new copy to the next lattice point
		copy->Shift(yshift + xshift + zshift);
		mols.push_back(copy);

		zshift += z_shift;
	  } 

	  xshift += x_shift;
	} 
	yshift += y_shift;
  }


  /* now that there are several copies of the molecules, let's go through and rename the atoms, and stick them all into a single molecule of sio2 */
  Molecule * qtz = new Molecule();
  qtz->Name("qtz");
  std::string new_name;
  int si=0; 
  int o=0; 
  int h=0;

  for (int mol = 0; mol < mols.size(); mol++)
  {
	for (int a = 0; a < mols[mol]->size(); a++)
	{
	  Atom * atom = mols[mol]->Atoms(a);
	  if (atom->Name().find("SI") != string::npos) {
		new_name = std::string("S") + itoa_base(si++, 36);
	  }
	  if (atom->Name().find("O") != string::npos) {
		new_name = std::string("O") + itoa_base(o++, 36);
	  }
	  if (atom->Name().find("H") != string::npos) {
		new_name = std::string("H") + itoa_base(h++, 36);
	  }
	  atom->Name(new_name);
	  qtz->AddAtom(atom);
	}
  }

  std::vector<Molecule *> final_mols;
  final_mols.push_back(qtz);
  PDBFile::WritePDB (final_mols);

return 0;
}


  /*
  // create the Beta-cstb molecule copy
  Molecule Beta-cstb_copy;
  Beta-cstb_copy.Name("Beta-cstb");
  RUN (Beta-cstb) {
  Beta-cstb_copy.AddAtom(Beta-cstb[i]);
  }

  // flip the Beta-cstb upside down
  VecR x2 (1,0,0);
  VecR y2 (0,1,0);
  VecR z2 (0,0,-1);

  Beta-cstb_copy.X(x2); Beta-cstb_copy.Y(y2); Beta-cstb_copy.Z(z2);
  MatR dcm = Beta-cstb_copy.DCMToLab(y);

  RUN (Beta-cstb_copy)
  Beta-cstb_copy[i]->Position (dcm * Beta-cstb_copy[i]->Position());

  // find oxygens above which to position a Beta-cstb
  vector <Atom *> Os;
  RUN (mols) {
  for (int j = 0; j < mols[i]->size(); j++) {
  Atom * o = (*mols[i])[j];
  if (o->Name() != "O") continue;
  if (o->Position()[y] < 3.63) continue;
  Os.push_back(o);
  }
  }

  // place a Beta-cstb above certain oxygens
  int iSecret, iGuess;
  // initialize random seed:
  srand ( time(NULL) );
  RUN (Os) {

// generate random number:
iSecret = rand() % 100 + 1;
if (iSecret > 90) continue;

// copy the starting molecule
Molecule * pm = new Molecule;
pm->Name("Beta-cstb");
RUN2 (Beta-cstb) {
Atom * pa = new Atom (*Beta-cstb_copy[j]);
pm->AddAtom(pa);
}

// shift the new Beta-cstb to where the oxygen is
Atom * si = pm->GetAtom("Si");
VecR shift = si->MinVector(Os[i]);
pm->Shift (shift + VecR(0,2.3,0));
mols.push_back(pm);
}
*/
