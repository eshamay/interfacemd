#include "mpisys.h"

#ifdef	SYSTEM_XYZ
MPIMolSystem::MPIMolSystem (int *argc, char ***argv, string xyzfile, VecR& syssize, string wannierpath) {
#endif
#ifdef	SYSTEM_AMBER
MPIMolSystem::MPIMolSystem (int *argc, char ***argv, string const prmtop, string const mdcrd, string const mdvel) {
#endif

	// First we set up the internal MPI data for each node
	_MPISystemInit(argc, argv);

	// And now set up all the actual working system of molecules
	#ifdef	SYSTEM_XYZ
	_MolSystemInit (xyzfile, syssize, wannierpath);
	#endif
	#ifdef	SYSTEM_AMBER
	_MolSystemInit (prmtop, mdcrd, mdvel);
	#endif

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
	delete[] _forces;

	if (_master) {
		delete _sys;
	}

return;
}

#ifdef	SYSTEM_XYZ
void MPIMolSystem::_MolSystemInit (string xyzfile, VecR& syssize, string wannierpath) {
#endif
#ifdef	SYSTEM_AMBER
void MPIMolSystem::_MolSystemInit (string prmtop, string mdcrd, string mdvel) {
#endif

	if (_master) {
		// first we establish the system on the master node for working the I/O
		#ifdef	SYSTEM_XYZ
		_sys = new XYZSystem (xyzfile, syssize, wannierpath);
		#endif
		#ifdef	SYSTEM_AMBER
		_sys = new AmberSystem (prmtop, mdcrd, mdvel);
		#endif
		_numAtoms = _sys->size();
		_numMols = _sys->Molecules().size();
	}
	// now we let each node know how many atoms we're working with

	MPI_Bcast (&_numAtoms, 1, MPI_INT, 0, _worldcomm);
	MPI_Bcast (&_numMols, 1, MPI_INT, 0, _worldcomm);

	// each atom will nead a place to store its position, velocity, and force data. Thus, every atoms will need 3 places (x,y,and z coords) for each vector.
	_positions 	= new double[_numAtoms * 3];
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
			const char * name = _sys->Molecules()[mol]->Name().c_str();
			while (*name!=(char)NULL) {
				molNames[molNamesOffset] = *name;
				name++;
				molNamesOffset++;
			}
			molNames[molNamesOffset] = (char)NULL;
			molNamesOffset++;
		}
	}

	// broadcast lotsa stuff for names of atoms and molecules
	MPI_Bcast (&atomNamesOffset, 1, MPI_INT, 0, _worldcomm);
	MPI_Bcast (&molNamesOffset, 1, MPI_INT, 0, _worldcomm);
	MPI_Bcast (atomNames, atomNamesOffset, MPI_CHAR, 0, _worldcomm);
	MPI_Bcast (molNames, molNamesOffset, MPI_CHAR, 0, _worldcomm);

	// then each node decomposes the names of the atoms and tacks them in to the system.
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
		_atoms[atom].SetMass();
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
			numAtoms[mol] = _sys->Molecules()[mol]->size();
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
	double pos[3];
	double force[3];
	VecR vpos, vforce;

	if (_master) {
		for (int atom = 0; atom < _numAtoms; atom++) {
			vpos = (*_sys)[atom]->Position();
			vforce = (*_sys)[atom]->Force();
			for (int i=0; i<3; i++) pos[i] = vpos[i];
			for (int i=0; i<3; i++) force[i] = vforce[i];
			//double * pos = (*_sys)[atom]->DPosition();
			//double * force = (*_sys)[atom]->DForce();
			_positions[atom*3] = pos[0];
			_positions[atom*3+1] = pos[1];
			_positions[atom*3+2] = pos[2];
			_forces[atom*3] = force[0];
			_forces[atom*3+1] = force[1];
			_forces[atom*3+2] = force[2];
		}
	}

	MPI_Barrier (_worldcomm);		// do we need to let everyone catch up here, or does Bcast act as a barrier?
	// bcast the positions and forces to everyone
	MPI_Bcast (_positions, _numAtoms*3, MPI_DOUBLE, 0, _worldcomm);
	MPI_Bcast (_forces, _numAtoms*3, MPI_DOUBLE, 0, _worldcomm);

	// now everyone sets their own positions and forces from the broadcast
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

// this should take a vector that was compiled on the master node, deconstruct it, and then ship it out to each node for reconstruction
void MPIMolSystem::BroadcastVector (vector<int>& vec) {

	int vec_size;

	// first we find out how big the vector is
	if (_master)		// note that this type of bcast is only from the master to other nodes
		vec_size = vec.size();

	MPI_Bcast (&vec_size, 1, MPI_INT, 0, _worldcomm);

	// then form the placeholder array on each node
	int ph[vec_size];

	// and the master populates the outgoing array
	if (_master) {
		RUN(vec)
			ph[i] = vec[i];
	}

	// master broadcasts to everyone else
	MPI_Bcast (ph, vec_size, MPI_INT, 0, _worldcomm);

	// and now the other nodes reconstruct the vector
	if (!_master) {
		for (int i = 0; i < vec_size; i++)
			vec[i] = ph[i];
	}

	// now each node should have the vector recreated from the master copy!

return;
}

void MPIMolSystem::UpdateHBonds (vector<Atom *>& int_atoms) {

/********** H-BONDS **********/
	// lastly, there are times when we need to know the H-bonded atoms to each atom in the system. Generally, this will only happen if the master is updating the connectivity matrix each time. Here, the master will go through and construct a list of all the Hbonded atoms (id's) for each atom in the system
	// the data structure is a series of atom ids, separated by a value of -1.
	std::vector<int> HBondIDs;

	// first the master (that has all the info) finds all the Hbonded atom IDs and packages them up
	if (_master) {
		RUN(int_atoms) {
			std::vector<Atom *> hbonds = int_atoms[i]->HBonds();

			RUN2 (hbonds) {
				HBondIDs.push_back (hbonds[j]->ID());
			}

			HBondIDs.push_back (-1);
		}
	}

	// now all the other nodes need to get in on this
	this->BroadcastVector (HBondIDs);


	vector<Atom *>::iterator atom;
	vector<int>::iterator hbond;
	atom = int_atoms.begin();
	hbond = HBondIDs.begin();

	if (!_master) {

		while (atom != int_atoms.end()) {
			// and then every node (sans the master) needs to clear out the HBond info from the last run
			(*atom)->ClearHBonds();

			// now take the new HBond data and put it all back together
			while (*hbond != -1) {
				(*atom)->FormHBond (_sys->Atoms(*hbond));
				hbond++;
			}
			hbond++;
			atom++;
		}
	}

return;
}
