#include "hno3.h"
int NitricAcid::numNitricAcids = 0;

NitricAcid::NitricAcid () : Molecule () {
	_centerofmass.Zero();
	_mass = 0.0;
	_name = "hno3";

	_n = (Atom *)NULL;
	_h = (Atom *)NULL;
	_oh = (Atom *)NULL;
	_o1 = (Atom *)NULL;
	_o2 = (Atom *)NULL;

	_set = false;

	++numNitricAcids;
}

NitricAcid::~NitricAcid () {
	--numNitricAcids;
}

/* To define the plane of a nitric acid molecule all we really need is the coordinates of the three oxygen atoms. Choosing one atom and finding the vectors from that atom to the other two, we can then find a vector normal to the plane of the molecule by finding the cross-product of the two vectors. That normal vector then defines the plane (sans a point of origin). */
VecR NitricAcid::MolecularPlaneVector () {

	// now we can calculate the two O-O vectors
	VecR oo1 = _o1->Position() - _oh->Position();
	VecR oo2 = _o2->Position() - _oh->Position();

	// lastly we find the normal vector to the molecular plane
	_molPlane = oo1 % oo2;

return(_molPlane.Unit());
}

// specialized routine for calculating the dipole of a nitric acid's NO2 group
bool NitricAcid::CalcNO2Dipole () {

	if(!_set) this->SetAtoms();
	this->UpdateCenterOfMass();

	// we've got all the atoms, now let's find the wannier centers that match those atoms of O1, O2, and N
	_no2wanniers.clear();
	for (VecR_it it = _wanniers.begin(); it != _wanniers.end(); it++) {

		// check for O1 wans
		double distance = (_o1->Position() - *it).Magnitude();
		if (distance < WANNIER_BOND) _no2wanniers.push_back (*it);

		// check for O2 wans
		distance = (_o2->Position() - *it).Magnitude();
		if (distance < WANNIER_BOND) _no2wanniers.push_back (*it);

		// check for N wans
		distance = (_n->Position() - *it).Magnitude();
		if (distance < WANNIER_BOND) _no2wanniers.push_back (*it);
	}

	// we've got all the atoms and wannier center locations... let's calculate the dipole
	_no2dipole.Zero();

	VecR temp = _centerofmass - _n->Position();
	_no2dipole += temp * 5.0;

	temp = _centerofmass - _o1->Position();
	_no2dipole += temp * 6.0;

	temp = _centerofmass - _o2->Position();
	_no2dipole += temp * 6.0;

	for (VecR_it it = _no2wanniers.begin(); it != _no2wanniers.end(); it++) {
		_no2dipole -= (_centerofmass - *it) * 2.0;
	}

return true;
}

void NitricAcid::PrintNO2 () const {

	cout << "Mag. " << _no2dipole.Magnitude() << "\tNum of Wanniers: " << _no2wanniers.size() << endl;
	_no2dipole.Print();
	_n->Print();
	_o1->Print();
	_o2->Print();

	for (VecR_it it = _no2wanniers.begin(); it != _no2wanniers.end(); it++) {
		printf ("wan) ");
		it->Print();
	}

return;
}

// The first time a nitric acid molecule is created, the atom pointers are set to point to the specific parts of the molecule
// this just helps speed up various calculations
void NitricAcid::SetAtoms () {
	if (!_set) {

		// here's the hydrogen and nitrogen atoms
		_h = (*this)["H"];
		_n = (*this)["N"];

		std::vector< std::vector<double> > oxygens;

		// we'll go through the three oxygens and find the one closest to the hydrogen - that becomes _oh. The other two then just get set as _o1 and _o2
		int i = 0;
		for (Atom_it it = _atoms.begin(); it != _atoms.end(); it++) {

			if ((*it)->Name().find("O") == string::npos) continue;

			std::vector<double> temp;
			// find the distance from each oxygen to the hydrogen
			double distance = (_h->Position() - (*it)->Position()).Magnitude();
			//double distance = _atoms[i]->Position().MinDistance (_h->Position(), Atom::Size());
			temp.push_back (distance);
			temp.push_back ((double)i++);

			oxygens.push_back (temp);
		}

		// sort them by distance
		sort(oxygens.begin(), oxygens.end());

		// and set the right oxygens
		_oh = _atoms[(int)(oxygens[0][1])];
		_o1 = _atoms[(int)(oxygens[1][1])];
		_o2 = _atoms[(int)(oxygens[2][1])];

		_set = true;
	}

	// lastly, let's set the OH vector (pointing from O to H)
	_voh = _h->Position() - _oh->Position();
	//_voh = _oh->Position().MinVector(_h->Position(), Atom::Size());

return;
}

VecR NitricAcid::NO2Bisector () {

	if (!_set) this->SetAtoms();

	VecR no1 = _o1->Position() - _n->Position();
	VecR no2 = _o2->Position() - _n->Position();

	VecR bisector = no1.Unit() + no2.Unit();

return (bisector);
}


int Nitrate::numNitrates = 0;

Nitrate::Nitrate () : Molecule () {
	_centerofmass.Zero();
	_mass = 0.0;
	_name = "no3";

	_n = (Atom *)NULL;
	_o1 = (Atom *)NULL;
	_o2 = (Atom *)NULL;
	_o3 = (Atom *)NULL;

	_set = false;

	++numNitrates;
}

Nitrate::~Nitrate () {
	--numNitrates;
}

// The first time a nitric acid molecule is created, the atom pointers are set to point to the specific parts of the molecule
// this just helps speed up various calculations
void Nitrate::SetAtoms () {

		// here's the hydrogen and nitrogen atoms
		_n = (*this)["N"];

		// now set the 3 oxygens
		std::vector<Atom *> oxygens;
		for (Atom_it it = _atoms.begin(); it != _atoms.end(); it++) {
			if ((*it)->Name() != "O") continue;
			oxygens.push_back (*it);
		}

		_o1 = oxygens[0];
		_o2 = oxygens[1];
		_o3 = oxygens[2];

		// while we're here we may as well also find the N-O bond vectors
		_no1 = _o1->Position() - _n->Position();
		_no2 = _o2->Position() - _n->Position();
		_no3 = _o3->Position() - _n->Position();
		//_no1 = _n->Position().MinVector(_o1->Position(), Atom::Size());

		_set = true;

return;
}

