#ifndef	MPISYS_H_
#define MPISYS_H_

// set this for dealing with the right type of system
#define SYSTEM_AMBER

#ifdef	SYSTEM_AMBER
#include "../ambersystem.h"
#endif
#ifdef	SYSTEM_XYZ
#include "../xyzsystem.h"
#endif

// because of the MPI-2 standard, we have to undef these (which are also included in the stdio.h)
#undef	SEEK_SET
#undef	SEEK_END
#undef	SEEK_CUR

#include "mpi.h"
#include <vector>

/*	This library provides functionality for running system analysis using mpi routines for tasks like SFG analysis, etc.
 *
 *	Be sure to define either:
 *		SYSTEM_XYZ
 *		or
 *		SYSTEM_AMBER
 *
 *	These choose between the two types of systems for analysis
 */
/* Define a few useful MPI macros */
#define BLOCK_LOW(id,p,n)	((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n)	(BLOCK_LOW((id)+1,(p),(n)) - 1)
#define BLOCK_SIZE(id,p,n)	(BLOCK_LOW((id)+1,(p),(n)) - BLOCK_LOW((id),(p),(n)))
#define BLOCK_OWNER(index,p,n)	(((p)*((index)+1)-1)/(n))

// this should hold all the data we want to work with to define our system of water molecules.
class MPIMolSystem {

private:
	#ifdef	SYSTEM_AMBER
	AmberSystem * _sys;			// the amber system object for the master node
	#endif
	#ifdef	SYSTEM_XYZ
	XYZSystem * _sys;			// The XYZ system, if we're working with one of those 
	#endif
	int			_numAtoms;			// number of atoms in the system - for all the nodes to know
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
	#ifdef SYSTEM_XYZ
	void _MolSystemInit (string xyzfile, VecR& syssize, string wannierpath);
	#endif
	#ifdef SYSTEM_AMBER
	void _MolSystemInit (string mdcrd, string mdcrd, string mdvel);
	#endif

	// Sets the molecule and atom names for the system
	void _InitAtoms ();

	// Initialize the molecules by adding the atoms into each molecule
	void _ParseMols ();

	// Update coordinates for position/velocity/forces, etc
	void _UpdateCoords ();

// Here's the meat and potatoes data that we have to work with
	double * _positions;
	double * _forces;

// Then we'll further abstract and form atoms and molecules to play with instead of just numbers for forces and positions
	std::vector<Atom> _atoms;
	std::vector<Molecule> _mols;

public:
	
	#ifdef	SYSTEM_XYZ
	MPIMolSystem (int *argc, char ***argv, string xyzfile, VecR& syssize, string wannierpath);
	#endif
	#ifdef	SYSTEM_AMBER
	MPIMolSystem (int *argc, char ***argv, string const prmtop, string const mdcrd, string const mdvel);
	#endif
	~MPIMolSystem ();

	bool Master () const { return _master; }
	int ID () const { return _id; }
	int Procs () const { return _p; }
	MPI_Comm WorldComm () const { return _worldcomm; }
	
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

	// system control commands
	void LoadNext ();
	#ifdef SYSTEM_AMBER
	void UpdateHBonds (vector<Atom *>& int_atoms);
	void UpdateGraph (vector<Atom *> int_atoms) { 
		_sys->UpdateGraph(int_atoms); 						// update the bonding data
		this->UpdateHBonds (int_atoms);						// and the Hbond data on each node
	}
	#endif
	
	std::vector<Molecule *>& Molecules () const { return _sys->Molecules(); }

	void BroadcastVector (vector<int>& vec);
};

#endif
