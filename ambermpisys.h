#ifndef AMBERMPISYS_H_
#define AMBERMPISYS_H_

#include "mpi.h"
#include "ambersystem.h"
#include <algorithm>

/* Define a few useful MPI macros */
#define BLOCK_LOW(id,p,n)	((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n)	(BLOCK_LOW((id)+1,(p),(n)) - 1)
#define BLOCK_SIZE(id,p,n)	(BLOCK_LOW((id)+1,(p),(n)) - BLOCK_LOW((id),(p),(n)))
#define BLOCK_OWNER(index,p,n)	(((p)*((index)+1)-1)/(n))

// this should hold all the data we want to work with to define our system of water molecules.
class MPIMolSystem {

private:
	AmberSystem * _sys;			// the amber system object for the master node
	int			_numAtoms;		// number of atoms in the system - for all the nodes to know
	int			_numMols;		// number of molecules in the system

// MPI-specific things
	MPI_Comm _worldcomm;
	MPI_Status _stat;

	int _id, _p;		// mpi rank of the process, and total number of processes
	bool _master;	// sets the alias for the master

	// functions for working with MPI
	void _MPISystemInit (int *argc, char ***argv);

	// final routines for shutting things down
	void _MPIFinalize ();

	// This will set up the molecule system with which we'll be working
	void _MolSystemInit (string prmtop, string mdcrd, string mdvel);

	// Sets the molecule and atom names for the system
	void _InitAtoms ();

	// Initialize the molecules by adding the atoms into each molecule
	void _ParseMols ();

	// Update coordinates for position/velocity/forces, etc
	void _UpdateCoords ();

// Here's the meat and potatoes data that we have to work with
	double * _positions;
	//double * _velocities;
	double * _forces;

// Then we'll further abstract and form atoms and molecules to play with instead of just numbers for forces and positions
	vector<Atom> _atoms;
	vector<Molecule> _mols;

public:
	
	MPIMolSystem (int *argc, char ***argv, string prmtop, string mdcrd, string mdvel);
	~MPIMolSystem ();

	bool Master () const 		{ return _master; }
	int ID () const 			{ return _id; }
	int Procs () const 			{ return _p; }
	MPI_Comm WorldComm () const { return _worldcomm; }
	MPI_Status * Stat () 		{ return &_stat; }

	
	// functions to divide up workloads
	int BlockLow (int n) { return BLOCK_LOW(_id, _p, n); }
	int BlockHigh (int n) { return BLOCK_HIGH(_id, _p, n); }
	int BlockSize (int n) { return BLOCK_SIZE(_id, _p, n); }

	// operator overloading
	Atom * Atoms (int index) { return &_atoms[index]; }
	Molecule * Molecules (int index) { return &_mols[index]; }

	// system info
	int NumMols () { return _numMols; }
	int NumAtoms () { return _numAtoms; }

	// Run through the system timesteps
	void LoadNext ();

	// find the location of the interfaces defined as the 50% points of the bulk number-densities.
	vector<double> FindInterfaces (string atomName, string residue);
};

MPIMolSystem::MPIMolSystem (int *argc, char ***argv, string prmtop, string mdcrd, string mdvel) {
	
	// First we set up the internal MPI data for each node
	_MPISystemInit(argc, argv);

	// And now set up all the actual working system of molecules
	_MolSystemInit (prmtop, mdcrd, mdvel);

	// set the atom and molecule names
	_InitAtoms ();

	// Build up the molecules by adding in the atoms
	_ParseMols ();

	// update the coords on the other nodes
	_UpdateCoords ();

}	

MPIMolSystem::~MPIMolSystem () {
	// before ending, let's do all the clean up we need!
	_MPIFinalize ();
}

// set up all things MPI
void MPIMolSystem::_MPISystemInit (int *argc, char ***argv) {
	// MPI definitions
	_worldcomm = MPI_COMM_WORLD;

	// Initialization of MPI
	MPI_Init(argc, argv);
	MPI_Comm_rank(_worldcomm, &_id);
	MPI_Comm_size(_worldcomm, &_p);

	_master = (_id == 0) ? true : false;		// set the master node

return;
}

void MPIMolSystem::_MPIFinalize () {

	// deallocation clean-up
	delete[] _positions;
	//delete[] _velocities;
	delete[] _forces;

	if (_master) {
		delete _sys;
	}

	MPI_Finalize();

return;
}

void MPIMolSystem::_MolSystemInit (string prmtop, string mdcrd, string mdvel) {

	if (_master) {
		// first we establish the system on the master node for working the I/O
		_sys = new AmberSystem (prmtop, mdcrd, mdvel);
		_numAtoms = _sys->size();
		_numMols = _sys->Molecules().size();
	}
	//
	// now we let each node know how many atoms we're working with
	MPI_Bcast (&_numAtoms, 1, MPI_INT, 0, _worldcomm);
	MPI_Bcast (&_numMols, 1, MPI_INT, 0, _worldcomm);

	// each atom will nead a place to store its position, velocity, and force data. Thus, every atoms will need 3 places (x,y,and z coords) for each vector.
	_positions 	= new double[_numAtoms * 3];
	//_velocities = new double[_numAtoms * 3];
	_forces 	= new double[_numAtoms * 3];

	_atoms.resize (_numAtoms);
	_mols.resize (_numMols);

return;
}

void MPIMolSystem::_InitAtoms () {

	// first thing we do is figure out the atom and molecule names - this part is easy
	//
	char atomNames[_numAtoms * 10];		// let's give ourselves some room to play with to set the names
	char molNames[_numMols * 10];
	int atomNamesOffset = 0, molNamesOffset = 0;

	// the master node gathers the names and serializes them for transmission
	if (_master) {
		for (int atom = 0; atom < _numAtoms; atom++) {
			const char * name = (*_sys)[atom]->Name().c_str();
			while (*name!=(char)NULL) {
				atomNames[atomNamesOffset] = *name;
				name++;
				atomNamesOffset++;
			}
			atomNames[atomNamesOffset] = (char)NULL;
			atomNamesOffset++;
		}
		for (int mol = 0; mol < _numMols; mol++) {
			const char * name = _sys->Molecules(mol)->Name().c_str();
			while (*name!=(char)NULL) {
				molNames[molNamesOffset] = *name;
				name++;
				molNamesOffset++;
			}
			molNames[molNamesOffset] = (char)NULL;
			molNamesOffset++;
		}	
	}

	// broadcast lotsa stuff for names...
	MPI_Bcast (&atomNamesOffset, 1, MPI_INT, 0, _worldcomm);
	MPI_Bcast (&molNamesOffset, 1, MPI_INT, 0, _worldcomm);
	MPI_Bcast (atomNames, atomNamesOffset, MPI_CHAR, 0, _worldcomm);
	MPI_Bcast (molNames, molNamesOffset, MPI_CHAR, 0, _worldcomm);

	// then each node decomposes the names of the atoms and tacks them in to the system
	// additionally, we'll do a few more things here (set the mass, id, etc)
	char * ptmp = atomNames;
	for (int atom = 0; atom < _numAtoms; atom++) {
		char name[10], *pname = name;
		for (int i = 0; i < 10; i++) name[i] = (char)NULL;

		// pull the name out of the list - unserialize it
		while (*ptmp!=(char)NULL) {
			*pname = *ptmp;
			pname++; ptmp++;
		}
		_atoms[atom].Name(string(name));	// and then each node sets each atom's name!
		ptmp++;

		// let's also set the mass
		_atoms[atom].SetMass(0.0);
		// and the atom's id, as well
		_atoms[atom].ID(atom);
		// let's also set the atom's position
		//_atoms[atom].Position().Set(_positions + atom * 3);
	}
	// now we do the molecules
	ptmp = molNames;
	for (int mol = 0; mol < _numMols; mol++) {
		char name[10], *pname = name;

		// pull the name out of the list - unserialize it
		while (*ptmp!=(char)NULL) {
			*pname = *ptmp;
			pname++; ptmp++;
		}
		_mols[mol].Name(string(name));	// and then set the molecule names :)
		ptmp++;

		// while we're at it, let's set the mol-ID, too
		_mols[mol].MolID(mol);
	}		

return;
}

void MPIMolSystem::_ParseMols () {

	// we need to set up the number of atoms in each molecule, so first the master node has to find that out
	int numAtoms[_numMols];
	if (_master) {
		for (int mol = 0; mol < _numMols; mol++) 
			numAtoms[mol] = _sys->Molecules(mol)->size();
	}

	// then the master let's all the others know about the breakdown of the atoms in the molecules
	MPI_Bcast (numAtoms, _numMols, MPI_INT, 0, _worldcomm);

	// then the others construct the molecules
	int offset = 0;
	for (int mol = 0; mol < _numMols; mol++) {
		for (int atom = 0; atom < numAtoms[mol]; atom++) {
			_mols[mol].AddAtom (&_atoms[offset]);
			offset++;
		}
	}
		
return;
}

void MPIMolSystem::_UpdateCoords () {
	
	// let's have the master grab all the updated positions and such
	if (_master) {
		for (int atom = 0; atom < _numAtoms; atom++) {
			double * pos = (*_sys)[atom]->DPosition();
			_positions[atom*3] = pos[0];
			_positions[atom*3+1] = pos[1];
			_positions[atom*3+2] = pos[2];

			double * force = (*_sys)[atom]->DForce();
			_forces[atom*3] = force[0];
			_forces[atom*3+1] = force[1];
			_forces[atom*3+2] = force[2];
		}
	}

	MPI_Bcast (_positions, _numAtoms*3, MPI_DOUBLE, 0, _worldcomm);
	MPI_Bcast (_forces, _numAtoms*3, MPI_DOUBLE, 0, _worldcomm);

	for (int atom = 0; atom < _numAtoms; atom++) {
		_atoms[atom].Position (_positions[atom*3], _positions[atom*3+1], _positions[atom*3+2]);
		_atoms[atom].Force (_forces[atom*3], _forces[atom*3+1], _forces[atom*3+2]);
	}

return;
}

// Routine to update all the positions, velocities, and forces
void MPIMolSystem::LoadNext () {

	if (_master) _sys->LoadNext();
	
	this->_UpdateCoords();

return;
}

vector<double> MPIMolSystem::FindInterfaces (string atomName, string residue) {

	vector<double> interfaces;

if (_master) {
	// the technique here will be to first construct the number-density function for the residue given. This is something that the master node can do on its own.
	
	double MIN =	-100.0;		// minimum of the histogram
	double MAX =	200.0;		// maximum of the histogram
	double	DZ =	0.10;		// bin size

	int numbins = (int)(MAX-MIN)/DZ;
	vector<int> histo (numbins, 0);

	for (int atom = 0; atom < _numAtoms; atom++) {
		
		// find the atom/residue
		if (!((*_sys)[atom]->Name() == atomName && (*_sys)[atom]->Residue() == residue)) continue;	

		int bin = (int) (((*_sys)[atom]->Y() - MIN) / DZ);

		histo[bin]++;
	}

	// now we have the histogram for the timestep. The next trick is to find the maxima. There are a couple scenarios:
	// 		1) The system lies within the boundaries of the box and so there are 2 interfaces
	// 		2) The system is split over the boundaries of the box, and so there will be 4 interfaces (virtually)
	// Finding the maximum, and then including all values within a certain statistical width (let's just say +/- 10%) should give us a good average value for the bulk. At that point we should be able to then find the 50% marks of the bulk values.
	
	double maxDensity = *max_element(histo.begin(), histo.end());

	double avgBulkDensity = 0.0;
	int num = 0;

	for (int i = 0; i < numbins; i++) {
		if (histo[i] > maxDensity * 0.75) {		// count only spots where the density is within 25% of the max
			num++;
			avgBulkDensity += histo[i];
		}
	}	

	avgBulkDensity /= num;

	cout << avgBulkDensity << endl;
	// now we have to find the closest position to the 50% mark on the number-density histogram.
	vector<double> bulk;
	for (int i = 0; i < numbins; i++) {
		if (histo[i] > avgBulkDensity * 0.5) 
			bulk.push_back (DZ*i+MIN);
	}
	sort(bulk.begin(), bulk.end());

	for (int i = 0; i < bulk.size(); i++) {
		if (i != bulk.size() - 1) {
			if (bulk[i+1] - bulk[i] > 10.0) {
				interfaces.push_back (bulk[i]);
				interfaces.push_back (bulk[i+1]);
			}
		}
	}
	if (bulk.size() == 0) {
		interfaces.push_back(bulk[0]);
		interfaces.push_back(bulk[bulk.size()-1]);
	}

for (int i = 0; i < interfaces.size(); i++) printf ("%f\t", interfaces[i]);
printf ("\n");

}
return interfaces;
}

#endif
