#ifndef WATERSFG_H_
#define WATERSFG_H_

/* This class performs various calculations in order to produce SFG spectra. The calculations all derive from the Morita and Hynes paper:

		A. Morita, J.T. Hynes, Chemical Physics 258 (2000) pp. 371-390

	The overall idea is that given a system trajectory for an aqueous surface, and the accompanying forces on each particle in the system, it is possible to produce SFG spectra for each timeframe. Various detailed calculations take place to find the shift in frequency of each OH-oscillator on each water molecule due to the forces it experiences from its environment. These shifts in frequency give rise to hyperpolarizabilities of the molecules, and the second-order susceptibility of the total (macroscopic) system. The intensity at each frequency is then calculated based on the chi-2 (second order susceptibility) magnitude. Review the paper above for further information on the procedure coded here.
*/


#include <complex>
#include <math.h>
#include "h2o.h"
#include "graph.h"

/***** Set this if using the dipole-dipole correction term *****/
#define DIPOLE_DIPOLE


/* system constants as defined in the morita-hynes paper */
/* Note: I think morita-hynes use atomic units (bohrs, hartrees, and the like) so we'll work in those*/
const double MPROT	=	1836.153;				// in atomic units, this is the mass of a proton. Mass of electron = 1.0
const double MOXY	=	(MPROT*8.0+8.0);		// mass of an oxygen atom in atomic units
const double MHYD	=	(MPROT*1.0+1.0);		// and mass of a hydrogen
const double M		=	(MOXY*MHYD)/(MOXY+MHYD);	// reduced mass of the OH bond/oscillator

//here's the constants given in the paper
const double k0	=	0.548;		// atomic units (force/length) Eh/ao/ao
const double l	=	-1.991;		// atomic units (force/length^2)	Eh/ao/ao/ao... yikes

// some very useful conversion factors
const double ANG2BOHR			=	1.889726125;					// angstroms to bohr radii
const double HARTREE2KCALPMOL	=	627.509;						// from hartree to kcal/mol
const double AMBER2ATOMIC		=	1.0/HARTREE2KCALPMOL/ANG2BOHR;	// convert amber forces (kcal/mol/A) into atomic force units

const double PREFACTOR	=	sqrt(k0/M)*(l/(2.0*k0*k0));		// the prefactor to multiply the bond force for freq shift (in atomic units) (eq 10c)
const double HZ2WAVENUMBER	=	3.335641e-11;				// convert from Hz to wavenumbers (cm-1)  (this is 1/c)
const double HZ2AU			=	2.418884324306202e-17;			// convert Hz to atomic units of frequency
const double AU2WAVENUMBER		=	HZ2WAVENUMBER/HZ2AU;			// convert from frequencies in atomic units to cm-1 (note: **not angular frequencies!** For that we need to fix the factor of 2*Pi)

const double UNCOUPLED_OH_FREQ	= 3706.5/AU2WAVENUMBER;			// the frequency of uncoupled OH bonds in the vapor phase (converted to frequency in atomic units)
const double COUPLING_CONST		= 49.5/AU2WAVENUMBER;				// Coupling const taken from the energy gap of the sym + antisym stretches (V12 in atomic units)

// Value of the magnitude of the dipole moment derivative (square root of the sum of the squares)
const double MU_DERIV_MAGNITUDE = sqrt(-0.058*-0.058 + 0.157*0.157);
// Length of a rigid SPC/E OH-bond
const double OH_LENGTH = 1.000*ANG2BOHR;	// in atomic units
const double OH_COM_LENGTH = MHYD*OH_LENGTH/(MHYD+MOXY);	// distance to the center of mass of the OH bond from the oxygen (atomic units of length)
const double MU_DERIV_LENGTH = MU_DERIV_MAGNITUDE * OH_COM_LENGTH;		// length of the differential dipole moment element used in calculating the dipole-dipole interaction energy (in atomic units)

/* Gamma, a "damping parameter" shows up as an arbitrary constant in the lorentzian part of the beta-spectrum. This number is tweaked to calibrate the spectra. DSW used a value of 2.0 and called it good after trying several values. */
const double GAMMA				= 2.0/AU2WAVENUMBER;				// the "damping parameter" first shown in Eq. 3 scaled down to some value (between 2 and 22 cm-1 as Dave put it)
const double GAMMA_SQ			= GAMMA*GAMMA;				// ...squared

// frequency limits for calculating the spectra (given in cm-1)
const double START_FREQ			= 2800.0/AU2WAVENUMBER;
const double END_FREQ			= 3800.0/AU2WAVENUMBER;
const double FREQ_STEP			= 1.0/AU2WAVENUMBER;		// step size when calculating the spectra
const double NUM_STEP			= (END_FREQ-START_FREQ)/FREQ_STEP;

class SFGCalculator {

private:

	//complex<double> I;

	bool 	_set;			// set when the water molecule has already gone through prelim calculations that don't need to be repeated for each calculation of Beta

	double _OH1FreqShift, _OH2FreqShift;		// The two frequency shift values calculated from the forces on the OH bonds
	double _w1, _w2;		// Shifted frequency values from the gas phase vibration frequency
	double _ws, _wa;		// The two normal mode (symmetric and antisymmetric) frequencies of the water molecule

	double _C1s, _C2s;		// The coefficients that are solved for in the matrix equation (Eq. 8) for the symmetric case
	double _C1a, _C2a;		// and the anti-symmetric case

	double _MuDerivA, _MuDerivS;			// the total mu derivative with the eigenvector weighting from Eq 7.
	double _AlphaDerivA, _AlphaDerivS;
	VecR MuDeriv1;				// dipole derivative vector from the paper
	VecR MuDeriv2;				// dipole derivative of the 2nd OH bond in the frame of the first
	MatR AlphaDeriv1;			// polarizability deriv matrix for the first OH bond in the frame of the 1st OH bond
	MatR AlphaDeriv2;			// polarizability deriv matrix of the 2nd OH in the frame of the 1st OH (through unitary transformation)

	std::vector< complex<double> > _Beta;		// hyperpolarizability of a given water (in the molecular frame)
	std::vector< complex<double> > _Chi;		// hyperpolarizability of a given water (in the molecular frame)

	MatR _DCM;				// rotation matrix for moving from the water molecular frame to the lab frame

	BondGraph * _graph;		// a connectivity matrix for analyzing water-bonding

	// calculate the dipole-dipole interaction potential between two dipoles (muA and muB) separated a distance R
	double DipolePotential (const VecR& muA, const VecR& muB, const VecR& R);

	double CouplingConstant (Water& water) const;

	void FreqShift (Water& water);		// Calculate the frequency shift on a given water-OH bond

	void WaterEigenSystem (Water& water);

	// calculate the total value of the collective product of the dipole and polarizability derivative terms
	void PolarizabilityAndDipoleDerivs (Water& water, const int s1, const int s2, const int p);
	//void PolarizabilityAndDipoleDerivs (Water& water);


public:

	//SFGCalculator (string polarization, coord axis);		// For loading up the entire system
	SFGCalculator (BondGraph * graph);		// For loading up the entire system

	void Reset () { _set = false; }
	// returns the summed beta rotated into the lab frame
	std::vector< std::complex<double> >& Beta (Water& water, const int s1, const int s2, const int p);
};

#endif
