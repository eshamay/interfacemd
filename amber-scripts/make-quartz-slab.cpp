#include "../pdbfile.h"
#include "quartz-utils.h"

using namespace std;

// Rename all the molecules into a single molecule
std::vector<Molecule *> UnitCellToSlab (Molecule *, VecR&, VecR&, VecR&);
Molecule * MakeSingleMolecule (std::vector<Molecule *>&);
void FixSingleMoleculeNaming (Molecule *);
void AddDanglingHydrogens (Molecule *);
void PrintMoleculePDB (Molecule *);

int main () {
  //PDBFile pdb ("quartz_100_short.pdb");
  PDBFile pdb ("beta-cristobalite.pdb");

  // create the quartz unitcell
  Molecule * uc = pdb.Molecules(0);

  //make copies of the quartz unitcell and shift them around
  //VecR x_shift (4.91, 0.0, 0.0);
  VecR x_shift (10.0, 0.0, 0.0);
  VecR y_shift (0.0, 10.0, 0.0);
  VecR z_shift (0.0, 0.0, 0.0);

  // Create the slab by creating and shifting copies of the unit cell
  Mol_ptr_vec mols = UnitCellToSlab(uc, x_shift, y_shift, z_shift);

  Molecule * mol = MakeSingleMolecule (mols);

  //AddDanglingHydrogens (mol);

  // Take care of renaming the atoms and the residue before printing it out
  FixSingleMoleculeNaming (mol);

  PrintMoleculePDB (mol);

  return 0;
}

// Makes a repeating slab out of a single unit cell by creating copies and shifting
Mol_ptr_vec UnitCellToSlab (Molecule * uc, VecR& x_shift,VecR& y_shift, VecR& z_shift) {

  std::vector<Molecule *> mols;
  VecR yshift, xshift, zshift;
  for (int i = 0; i < 1; i++) {
    zshift.Zero();
    yshift.Zero();

    for (int j = 0; j < 6; j++) {
      xshift.Zero();

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

	xshift += x_shift;
      } 

      yshift += y_shift;
    } 
    zshift += z_shift;
  }

  return mols;
}

/* Take a set of molecules and return a single molecule with all the atoms */
Molecule * MakeSingleMolecule (Mol_ptr_vec& mols) {
  Molecule * mol = new Molecule();

  for (Mol_it mols_i = mols.begin(); mols_i != mols.end(); mols_i++) {

    for (Atom_it atom_i = (*mols_i)->begin(); atom_i != (*mols_i)->end(); atom_i++) {

      mol->AddAtom((*atom_i));
    }
  }

  return mol;
}

/* Single molecule naming - takes a group of residues and sets the naming scheme so that the output is a single residue with lots of renamed atoms */
void FixSingleMoleculeNaming (Molecule * mol) {
  /* now that there are several copies of the molecules, let's go through and rename the atoms, and stick them all into a single molecule of sio2 */
  mol->Name("qtz");
  std::string new_name;
  int si=0; 
  int o=0; 
  int h=0;

  for (Atom_it atom_i = mol->begin(); atom_i != mol->end(); atom_i++) {

    if ((*atom_i)->Name().find("SI") != string::npos) {
      new_name = std::string("S") + itoa_base(si++, 36);
    }
    if ((*atom_i)->Name().find("O") != string::npos) {
      new_name = std::string("O") + itoa_base(o++, 36);
    }
    if ((*atom_i)->Name().find("H") != string::npos) {
      new_name = std::string("H") + itoa_base(h++, 36);
    }

    (*atom_i)->Name(new_name);
  }

  return;
}

void AddDanglingHydrogens (Molecule * mol) {

  Atom_ptr_vec atoms = mol->Atoms();
  Atom_ptr_vec Hs;

  /* cycle through each atom in order of the pdb file */
  for (Atom_it atom_i = mol->begin(); atom_i != mol->end(); atom_i++) {
    Atom * atom = *atom_i;
    VecR pos = atom->Position();

    /* find the surface Oxygens */
    if (atom->Name().find("O") != string::npos && pos[z] > -8.5)
    {
      /* add H's to the surface */
      Atom * h = new Atom ("H", pos + VecR (0.0, 0.0, 1.1));
      Hs.push_back(h);
    }
    /* oxygens on the slab's side */
    if (atom->Name().find("O") != string::npos && pos[y] > 36.0) {
      Atom * h = new Atom ("H", pos + VecR (0.0, 1.1, 0.0));
      Hs.push_back(h);
    }
  }

  // Add the Hs into the molecule
  RUN (Hs)
    mol->AddAtom(Hs[i]);

  return;
}

void PrintMoleculePDB (Molecule * mol) {
  Mol_ptr_vec final_mols;
  final_mols.push_back(mol);
  PDBFile::WritePDB (final_mols);
  return;
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
