//#include "../watersystem.h"
#include "../xyzsystem.h"
#include "../pdbfile.h"

using namespace std;

int main () {
	XYZFile xyz ("decane.xyz");
	Atom::Size (VecR(40.0, 40.0, 40.0));

	Molecule mol;
	mol.Name("dec");
	mol.MolID(0);

	RUN (xyz)
		mol.AddAtom (xyz[i]);

	Atom * c1 = mol.GetAtom("C1");
	Atom * c10 = mol.GetAtom("C10");

	// the molecular axis from the C1 to the furthest C10
	VecR _z = c10->MinVector(c1).Unit();

	// To find the x axis - take a hydrogen off the molecule, and find the minvector from the Si to the H. Then take out any z-component, and grab the unit vector.
	Atom * h1 = mol.GetAtom("H30");

	VecR _x = c10->MinVector(h1);
	_x = _x - (_z * (_x * _z));
	_x = _x.Unit();

	// y is perpendicular to x and z...
	VecR _y = (_z % _x).Unit();

	mol.X(_x);
	mol.Y(_y);
	mol.Z(_z);

	MatR dcm = mol.DCMToLab(y);

	VecR X (1, 0, 0);
	VecR Y (0, 1, 0);
	VecR Z (0, 0, 1);

	Atom * pa;
	RUN(mol) {
		pa = mol[i];
		pa->Position(dcm * pa->Position());
	}

	//mol.Print();
	std::vector<Molecule *> mols;
	mols.push_back(&mol);
	PDBFile::WritePDB (mols);

return 0;
}
