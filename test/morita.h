#ifndef MORITA_H_
#define MORITA_H_

#include <stdio.h>
#include <stdlib.h>
#include "../utility.h"
#include "../watersfg.h"
#include "../matrixr.h"

// if we're running this on an MPI system
//#define MPI_SYS

#ifdef	MPI_SYS
#include "../mpi/mpisys.h"
#else
#include "../ambersystem.h"
#endif

#define PRMTOP	"prmtop"
#define MDCRD	"mdcrd"
#define FORCE	"mdvel"


// set the interface correctly - different simulation = different location!!
#define INTERFACE_LOW		68.0			// the location of the interface (the top one, only, for now)
#define INTERFACE_HIGH		90.0
#define PBC_FLIP			15.0			// used for funcy periodic boundaries

const coord axis = y;
#define	 OUTPUT_FREQ	10					// how often the output file will be written (# of timesteps/10)
#define	 TIMESTEPS		20000				// # of timesteps to process through
#define	 POLARIZATION	"SPS"				// the polarization scheme used for the analysis
#define	 OUTPUTFILE		"morita.sfg.dat"	// name of the output file for the final spectra


#ifdef MPI_SYS
void OutputHeader (MPIMolSystem& sys);
void OutputStatus (int const count, MPIMolSystem& sys);
void MPI_PackageChi (std::vector< std::complex<double> >& TimestepChi, std::vector< std::complex<double> >& TotalChi, MPIMolSystem& sys);
void FindInterfacialAtoms (vector<Molecule *>& int_mols, vector<Atom *>& int_atoms, MPIMolSystem& sys);
#else
void OutputHeader ();
void OutputStatus (int const count);
void FindInterfacialAtoms (vector<Water *>& int_mols, vector<Atom *>& int_atoms, AmberSystem& sys);
#endif

void OutputData (FILE * fp, vector< complex<double> >& chi);

void CollectChi (std::vector< std::complex<double> >& MolChi, std::vector< std::complex<double> >& TotalChi);


#endif
