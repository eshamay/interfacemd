#define OH_BOND	1.3
#define H_BOND	2.46 
#define WANNIER_BOND	0.6			// max distance that a wannier center will sit away from an Oxygen on a water molecule

#include <vector>
#include "vecr.h"
#include "xyzfile.h"
#include "wannier.h"
#include "molecule.h"
#include "utility.h"
	

// For naming the bond-types, let's initialize the COORD_NAMES array;
string COORD_NAMES[12] = {"Unbound", "O", "OH", "OHH", "OO", "OOH", "OOHH", "H", "HH", "OOOHH", "OOHHH", "OVER"};

// Encoding of the different coordination types
enum coordination {UNBOUND, O, OH, OHH, OO, OOH, OOHH, H, HH, OOOHH, OOHHH, OVER};

// A non-member function for getting a nice string representation of a coordination
const string CoordName (coordination coord)	{ return COORD_NAMES[coord]; }

/* note that the functions below are geared for systems with water. The dipoles and coordinations will run through the system, find the 
water molecules, and then calculate as necessary for the waters. All other molecules are left untouched. The output of the functions
is generally a vector that has one entry for each water. The entries are based on the Oxygens as there is only one per water. Thus,
the vector _oxygens has a size corresponding to the number of waters in the system, and each entry is the index of the O-atom in the
atom list _atoms
*/
class WaterSystem {
	
	XYZFile		 			_atoms;			// Atomlist parsed from an xyz file
	vector<int>				_oxygens;		// holds the indexes of the water oxygens
	vector<VecR>			_dipoles;		// An array that holds the dipoles of the waters from wannier center calcs
	vector< vector<double> > 	_waterMap;	// A connectivity matrix that holds water connection data for efficient calcs
	vector< vector<int> >	_hBonds;		// Storage for the number of h-bonds each atom makes (includes both O and H atoms)
											/* This vector has a length of the number of atoms in the system, and 2 elements each to hold
											   the number of acceptor-bonds, and the number of donor-bonds.
											   This, the vector is constructed as: _hBonds[atom][0] is the #acceptor-bonds,
											   									   _hBonds[atom][1] is the #donor-bonds.	*/
	vector<coordination>	_coords;		// Gives the coordination for each water at each timestep
	vector<VecR>			_bisectors;		// A list of all the bisector vectors of water
	vector<double>			_rotations;		// A list of h2o-rotation angles in each molecule
	WannierCenters			_wanniers;		// the wannier centers of the water system (not ordered)
	VecR					_com;			// center of mass of the system

	bool _mapUpdated;						// Has the map been updated for the current timestep?
	bool _wanniersloaded;

	VecR				_size;				// A vector holding the system size

	void _FindHBonds ();					// Find the number of H-bonds each atom is involved in

public:
	// constructors
	WaterSystem (string filepath, VecR size, string wannierpath);
	
	// Controller & Calculation methods
	void UpdateConnectivityMatrix ();
	void FindCoordinations ();				// Determine the H-bond coordination type for each water in the system
	void FindBisectorsAndRotations ();		// go through each water and find it's bisector vector
	void FindDipoles ();					// Calculate the dipole vector for every water
	void FindCenterOfMass ();

	void LoadNext ();	 					// Update the system to the next timestep 
	void LoadFirst ();
	void LoadLast ();
	void Seek (int step);

	// Input
	void Size (VecR size) { _size = size; }		// Set the system size
	VecR Size () const { return _size; }		// returns the system size. Not the uppercase S to distinguish from the number of atoms

	// Output
	const XYZFile& 	AtomList () 		{ return _atoms; 			}
	int 		 	size ()				const { return _atoms.size();		}

	const int First ()					const { return _atoms.First();	}
	const int Current ()				const { return _atoms.Current(); 	}
	const int Last ()					const { return _atoms.Last(); 	}

	vector<coordination>&	Coordinations (){ return _coords; 			}		// just return the current _coords
	vector<VecR>& 			Bisectors () 	{ return _bisectors;		}		// give back the listing of vector bisectors
	vector<double>& 		Rotations () 	{ return _rotations;		}
	vector< vector<double> >& 	Map ()		{ return _waterMap;			}
	vector<VecR>& 			Dipoles () 		{ return _dipoles;			}
	vector<int>& 			Oxygens ()		{ return _oxygens;			}
	VecR& 					CenterOfMass ()	{ return _com;				}
	WannierCenters&			Wanniers ()		{ return _wanniers;			}

	// operators
	Atom& operator[] (int index) { return _atoms[index]; }
};

WaterSystem::WaterSystem (string filepath, VecR size, string wannierpath = "") {
	
	_atoms = XYZFile (filepath);			// initialize the file to get things rolling
	
	if (wannierpath != "") {
		_wanniers = WannierCenters (wannierpath);
		_wanniersloaded = true;
	}
	else _wanniersloaded = false;

	// set the system size
	_size = size;
	Atom::Size (size);

	// Initialize the waterMap (2-d connectivity matrix for atom-atom bonds)
	_waterMap = vector< vector<double> > (_atoms.size(), vector<double> (_atoms.size(), 0.0));

	// then update the connectivity matrix
	_mapUpdated = false;
	this->UpdateConnectivityMatrix();

	_hBonds = vector< vector<int> > (_atoms.size(), vector<int> (2, 0));

	_bisectors = vector<VecR> ();

	_rotations = vector<double> ();
}

void WaterSystem::UpdateConnectivityMatrix() {

	double distance;

    /* map elements store the connections between oxygens and hydrogens in a water system. If a map value is 0.0, then there is no connection between the two (meaning that they are neither hydrogen bonded or OH bound).
		If an oxygen-hydrogen pair are within a type of bonding distance, then one of two things can happen.
			a) for OH-bonds, the actual value of the bond-length is stored in the map element.
			b) for H-bonds, the NEGATIVE of the bond length (distance from the O to the H) is stored.
			c) as explained above, a value of 0.0 is used if the two atoms are not bonded (OH or H-bond)
	*/

	// clear out the map to start out
	RUN (_waterMap) {
		RUN2 (_waterMap) {
			_waterMap[i][j] = 0.0;
		}
	}

	// run through each particle-pair 
	RUN (_waterMap) {
		RUN2 (_waterMap) {

			// now let's see if we have any bonding

			if (_atoms[i].Name() != "O" || _atoms[j].Name() != "H") continue;
			// Note: we're only checking here for bonds between O's and H's... for H-bonding, other elements may be involved!
			// this may need to be changed at some point.

			// let's get the distance between the two atoms
			distance = _atoms[i] - _atoms[j];

			// determine if we have an OH-bond (compare to the OH_BOND definition) and mark those as -1
			if (distance < OH_BOND)	{
				_waterMap[i][j] = distance;
				_waterMap[j][i] = distance;
				_waterMap[i][i] += 1.0;
			}

			// now let's see if we have an H-bond, and mark it as -2 if it falls in the bond-length range
			else if (distance > OH_BOND && distance <= H_BOND) {
				_waterMap[i][j] = -distance;
				_waterMap[j][i] = -distance;
			}
		}
	}

	_mapUpdated = true;
}

void WaterSystem::_FindHBonds () {
	
	// first zero out the _hBond array
	RUN (_hBonds) {
		_hBonds[i][0] = 0;
		_hBonds[i][1] = 0;
	}

	// Run through each particle pair, and count the H-bonds
	RUN (_waterMap) {
			
		if (_atoms[i].Name() != "O" || _waterMap[i][i] != 2.00000) continue; // only look at O atoms that are in proper waters (H2O)

		/* There are two things to do for each oxygen in the watermap. We will first take a trip over the entire row of the map and find every bound hydrogen. Once we find the index of the hydrogen we'll take a very brief detour to the hydrogen's row of the map to count up the number of donor-bonds it forms. Then we return to the oxygen row and count its own acceptor bonds formed. All these bond numbers will be stored in the _hBond vector so that we may easily figure out the coordination of the molecule. */

		// for each O, let's run through the entire row of the watermap and find some things
		RUN2 (_waterMap[i]) {
			if (i == j) continue;		// skip diagonal elements of the map

			// next, let's count the number of acceptor bonds the oxygen forms
			if (_waterMap[i][j] < 0.0) {
				_hBonds[i][0]++;
			}

			// let's find the two OH-bound hydrogens and count the number of donor bonds they form
			if (_waterMap[i][j] > 0.0) {
				// if a bound-hydrogen is found, then we briefly leave the oxygen row to sojourn over the hydrogen's row and count up the number of H-bonds it makes. 
				for (int k = 0; k < _waterMap[j].size(); k++) {
					if (_waterMap[j][k] < 0.0) _hBonds[i][1]++;		// the water's h-bonding info is stored according to the oxygen's index
				}
			}

		}
	}

	// and that's that! we've now got a complete listing of each water-molecule's acceptor and donor bond numbers
}

/* Here we construct a vector (1D) with an entry for every atom in the system. For each of the O's we check to see if it's
   a water O (by counting the number of OH bonds - diagonal elements of the connectivity matrix). For each water, we count
   the number of acceptor and donor bonds from the O's and H's, and then classify it according to coordination type.
			now there exist 9 possible coordinations for an H-bonded water:
				1: no bonds
				2: O
				3: OH
				4: OHH
				5: OO
				6: OOH
				7: OOHH
				8: H
				9: HH
				10: OOOHH	(note: both of these are over-coordinated waters)
				11: OOHHH
*/
void WaterSystem::FindCoordinations () {
	
	// first refresh the h-bond listing
	if (!_mapUpdated) this->UpdateConnectivityMatrix();
	this->_FindHBonds();

	// and then clear out the previous listing
	_coords.clear();
	_oxygens.clear();

	// going through each atom in the system
	RUN (_waterMap) {

		if (_atoms[i].Name() != "O" || _waterMap[i][i] != 2.00000) continue;		// we'll only look at the proper waters for now

		_oxygens.push_back(i);		// here we update the _oxygens to keep track of where our 'proper' waters are

		// now comes the ugly part of assigning each bonding type a name (as outlined above)

		// 1: Uncoordinated
		if (_hBonds[i][0] + _hBonds[i][1] == 0) { 
			_coords.push_back(UNBOUND);
			continue;
		}

		//2: O
		if (_hBonds[i][0] == 1 && _hBonds[i][1] == 0) {
			_coords.push_back(O);
			continue;
		}

		//3: OH
		if (_hBonds[i][0] == 1 && _hBonds[i][1] == 1) {
			_coords.push_back(OH);
			continue;
		}

		//4: OHH
		if (_hBonds[i][0] == 1 && _hBonds[i][1] == 2) {
			_coords.push_back(OHH);
			continue;
		}

		//5: OO
		if (_hBonds[i][0] == 2 && _hBonds[i][1] == 0) {
			_coords.push_back(OO);
			continue;
		}

		//6: OOH
		if (_hBonds[i][0] == 2 && _hBonds[i][1] == 1) {
			_coords.push_back(OOH);
			continue;
		}

		//7: OOHH
		if (_hBonds[i][0] == 2 && _hBonds[i][1] == 2) {
			_coords.push_back(OOHH);
			continue;
		}

		//8: H
		if (_hBonds[i][0] == 0 && _hBonds[i][1] == 1) {
			_coords.push_back(H);
			continue;
		}

		//9: HH
		if (_hBonds[i][0] == 0 && _hBonds[i][1] == 2) {
			_coords.push_back(HH);
			continue;
		}

		//10: OOOHH - over-coordinated species
		if (_hBonds[i][0] == 3 && _hBonds[i][1] == 2) {
			_coords.push_back(OOOHH);
			continue;
		}

		//11: OOHHH - over-coordinated species
		if (_hBonds[i][0] == 2 && _hBonds[i][1] == 3) {
			_coords.push_back(OOHHH);
			continue;
		}

		//12: 6-coordinated species!!!
		if (_hBonds[i][0] + _hBonds[i][1] > 5) {
			_coords.push_back(OVER);
			continue;
		}
	}
}

void WaterSystem::FindBisectorsAndRotations () {

	/* First we need to run through to each Oxygen and see if it's an h2o oxygen or not. If it is, then we grab both the 
	   OH vectors and normalize them, then add the two norm(OH) vectors to get the bisector.
	*/
	VecR z (0.0, 0.0, 1.0);		// the z-axis

	_bisectors.clear();
	_rotations.clear();
	_oxygens.clear();
	if (!_mapUpdated) this->UpdateConnectivityMatrix();

	// look through the watermap at each oxygen and find the number of bound hydrogens
	RUN (_atoms) {
		int h1 = 0, h2 = 0;		// indices for the bound hydrogens
		
		// look for O's with only 2 bound H's (for now)
		if (_atoms[i].Name() == "O" && _waterMap[i][i] == 2.00000) {
			
			_oxygens.push_back(i);

			for (int j = 0; j < _waterMap.size(); j++) {
				// find the bound H's
				if (_waterMap[i][j] > 0.0 && i != j) {
					// load them into the two indices
					if (!h1) h1 = j;
					else if (!h2) h2 = j;
				}
			}

			// now calculate the bisector vector (note that both vectors are normalized (.Unit())
			VecR vh1 = _atoms[h1].Position();
			VecR vh2 = _atoms[h2].Position();
			VecR vO = _atoms[i].Position();
			_bisectors.push_back( (vh1 - vO).Unit() + (vh2 - vO).Unit() );

/* There are two vectors that are needed to describe a water molecule's orientation. The first is the bisector which can give information about the molecule's alignment with the z-axis. However, although the aligned in any given direction (with an angle ranging from 0 (aligned with the z-axis) to pi/2 (perpendicular) to pi (aligned in the negative z-direction), one more angle is needed to describe the orientation with respect to rotation about the bisector. To note: it is possible for a water molecule to lie with bisector parallel to the surface of the water and both OH bonds lying in the plane, or both tilted out of plane (i.e. one pointing into the bulk, and one out of the surface).

The second vector needed is one describing the plane of the molecule and it's relation to the Z-axis. To define this angle, we first find the vector normal to the bisector and the z-axis. Then we take the dot product of the plane-vector with the molecule's H-H vector to get a value between 0 and pi for the molecule's rotation in the surface-normal plane.
*/
			// now calculate the bisector-z-axis normal vector (molecular plane-vector)
			VecR plane = _bisectors[i] % z; // (cross-product of the bisector and z-axis)
			
			VecR vhh = vh1 - vh2;

			// Then calculate the (cos)angle made between plane vector and the HH vector;
			_rotations.push_back ( ( vhh * plane ) * (1.0/vhh.Magnitude()) * (1.0/plane.Magnitude()) );
		}
	}
}

void WaterSystem::FindDipoles () {

	// we need a wannier-center file in order to process dipole information!
	if (!_wanniersloaded) {
		printf ("\nTrying to calculate dipoles without a proper wannier-center file!\n\n");
		exit(1);
	}

	/* a few steps to calculating the dipole for each water molecule. First we find the molecule with O, and both H's.
	   Then we run through the wannier centers and find the 4 closest centers to each oxygen. The dipole is then a function
	   of those vectors as follows:
	   	vDipole = 6*vO + vH1 + vH2 - 2*sum(vWannier)
	   Straightforward enough. The vectors for the dipoles will be stored in the list _dipoles. It has to be updated 
	   each timeframe where we want the dipoles.

	   // so first let's rock and roll on finding the vectors for the O, h1, h2, and the centers
	 */

	_dipoles.clear();
	_oxygens.clear();
	if (!_mapUpdated) this->UpdateConnectivityMatrix();

	RUN (_atoms) {
		if (_atoms[i].Name() != "O" || _waterMap[i][i] != 2.00000) continue;	// only check proper waters

		_oxygens.push_back(i);
		vector<VecR> wanniers = _wanniers.Centers();

		VecR vOxy = _atoms[i].Position();
		vector<VecR>	hyds;	// hydrogens of the water
		vector<VecR>	wans;	// associated wannier centers
		vector<int>		wanind;

		// grab the hydrogens from the atom list
		RUN2 (_waterMap) {
			if (_waterMap[i][j] > 0.0 && i != j) {
				_atoms[j].Wrap(vOxy);
				hyds.push_back (_atoms[j].Position());
			}
		}

		// then load up all the nearby wannier centers
		double distance;
		RUN2 (_wanniers) {
			distance = vOxy.MinDistance (wanniers[j], _size);

			// Check to see if the center is near the water
			if (distance < WANNIER_BOND) {
				wanniers[j].Wrap(_size, vOxy);
				wans.push_back (wanniers[j]);
				wanind.push_back(j);
			}
		}

		// Okay, we have all the centers and h's accounted for! Let's do some math to get the dipoles
		// calculate the dipole based on the position vectors of the atoms and wannier centers 
		//dipole = ( vO * 6.0 ) + hyds[0] + hyds[1] - ( (wanniers[0] + wanniers[1] + wanniers[2] + wanniers[3])*2.0);
	
		if (hyds.size() != 2 || wans.size() != 4) {
			printf ("[%d  %d  %d] Oxygen::%d\n", this->Current(), hyds.size(), wans.size(), i);

			RUN2(hyds) {
				printf ("H ::\t\t");
				hyds[j].PrintVector();
			}
			RUN2(wans) {
				printf ("W(%d) ::\t\t", wanind[j]);
				_wanniers[wanind[j]].PrintVector();
			}
		}
		else {
			VecR dipole;

			dipole = (vOxy * 6.0) + hyds[0] + hyds[1] - ( (wans[0] + wans[1] + wans[2] + wans[3])*2.0);

			if (dipole.Magnitude() > 3.0) {
				printf ("magnitude (% .4f) [%d  %d]\n", dipole.Magnitude(), this->Current(), i);
			}
			else _dipoles.push_back (dipole);	// and add the dipole to the list
		}
	}
}

void WaterSystem::LoadFirst () {
	_mapUpdated = false;
	_atoms.LoadFirst();
	
	if (_wanniersloaded) _wanniers.LoadFirst();
}

void WaterSystem::LoadLast () {
	_mapUpdated = false;
	_atoms.LoadLast();
	if (_wanniersloaded) _wanniers.LoadLast();
}

void WaterSystem::Seek (int step) {
	_mapUpdated = false;
	_atoms.Seek (step);
	if (_wanniersloaded) _wanniers.Seek (step);
}

void WaterSystem::LoadNext () {
	_mapUpdated = false;
	_atoms.LoadNext();
	if (_wanniersloaded) _wanniers.LoadNext();
}
