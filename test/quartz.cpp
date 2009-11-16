//#include "../watersystem.h"
#include "../xyzsystem.h"
#include "../pdbfile.h"

using namespace std;

int main () {
  XYZFile xyz ("beta-cristobalite-unitcell.xyz");
  //  XYZFile pds ("beta-cristabolite-slab.xyz");
  Atom::Size (VecR(90.0, 90.0, 90.0));

  // create the quartz unitcell
  Molecule uc;
  uc.Name("qtz");
  RUN (xyz) {
    uc.AddAtom(xyz[i]);
  }

  //make copies of the quartz unitcell and shift them around
  VecR x_shift (0.0, 3.583, 3.583);
  VecR y_shift (3.583, 0.0, 3.583);
  VecR z_shift (3.583, 3.583, 0.0);

  vector<Molecule *> mols;
  int count = 9;
  for (int i = 0; i < count; i++) {
    VecR shift;
    shift += z_shift * (double)i;
    for (int j = 0; j < count; j++) {
      shift += x_shift;    
      for (int k = 0; k < count; k++) {
	shift += y_shift;

	// make a copy of the unitcell
	Molecule * copy = new Molecule;
	copy->Name("qtz");
	// make copies of all the atoms in the unitcell
	for (int atom = 0; atom < uc.size(); atom++) {
	  Atom * pa = new Atom (*uc[atom]);
	  copy->AddAtom(pa);
	}

	// shift the new copy to the next lattice point
	copy->Shift(shift);
	mols.push_back(copy);
      } 
      shift -= y_shift * (double)count;
    } 
    shift -= x_shift * (double)count;
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
  PDBFile::WritePDB (mols);
  return 0;
}
