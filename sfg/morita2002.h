#ifndef MORITA2002_H_
#define MORITA2002_H_

#include "../analysis.h"
#include "sfgunits.h"
#include "../tensor.h"
#include <Eigen/LU>
//#include <boost/shared_ptr.hpp>
#include <boost/timer.hpp>
#include <iostream>
//#include <fftw3.h>
#include <mkl_blas.h>
#include <mkl_lapack.h>


namespace morita {

  USING_PART_OF_NAMESPACE_EIGEN

  /*	From the morita 2002 water model
		const double BONDLENGTH_EQ	= 0.9575;		// in angstroms
		const double ANGLE_EQ 		= 104.51*M_PI/180.0;		// in radians
		const double CHARGE_H_EQ		= 0.3285;		// charge units (atomic units?)
		const double CHARGE_O_EQ		= -0.6570;
   */

  /* from the equilibrium values of the spcfw model */
  const double BONDLENGTH_EQ	= 1.012 * sfg_units::ANG2BOHR;		// in atomic units
  const double ANGLE_EQ 		= 113.24;		// in degress
  const double CHARGE_H_EQ		= 0.41;		// charge units (atomic units?)
  const double CHARGE_O_EQ		= -0.82;


  // these are the constants from the 2002 Morita work for calculating the dipole moment
  const double C1 =  0.1396,
		C2 = -0.1196,
		C3 = -0.0164,
		C4 = -0.0288,
		C5 = -0.0516,
		C6 =  0.0532,
		C7 = -0.0699,
		C8 =  0.0169,
		C9 =  0.1142;

  const double D1 = 4.6077,
		D2 = 4.8894,
		D3 = 5.5062,
		D4 = 1.6890,
		D5 = 1.6102,
		D6 = 7.3812,
		D7 = 3.4710;



  // a new water that has all the needed pieces for our calculations
  class MoritaH2O : public Water {
	public:

	  MoritaH2O (const Molecule& molecule) 
		: Water(molecule) { return; }
	  MoritaH2O (const Molecule * molecule) 
		: Water(*molecule) { return; }

	  void SetDipoleMoment ();
	  void SetPolarizability ();

	  MatR& Polarizability() { return _alpha; }

	protected:

	  MatR _alpha1, _alpha2;	// polarizabilities of the two OH bonds
	  MatR _alpha;

	  double dR1, dR2, dA;		// displacements of the OH bonds and the HOH angle from their equilibrium values
	  double X1, X2;
	  double qO, qH1, qH2;		// charges calculated from the MoritaHynes 2002 method
  };	// morita-h2o

  typedef MoritaH2O * MoritaH2O_ptr;
  typedef std::vector<MoritaH2O_ptr> Morita_ptr_vec;
  typedef Morita_ptr_vec::const_iterator Morita_it;




  // dipole field tensor 'T' as used in the morita&hynes paper
  // it's a square 3Nx3N matrix, where N = number of particles
  class DipoleFieldTensor : public MatR {
	public:
	  DipoleFieldTensor (const MoritaH2O_ptr wat1, const MoritaH2O_ptr wat2);
  }; // Dipole field tensor



  class SFGAnalyzer : public Analyzer<AmberSystem> {
	public:
	  SFGAnalyzer (WaterSystemParams& params);
	  ~SFGAnalyzer ();

	  void Setup ();
	  void Analysis ();
	  void DataOutput ();
	  void PostAnalysis ();

	private:
	  Morita_ptr_vec		all_wats;
	  Morita_ptr_vec		cutoff_wats;
	  VectorXd				_p;
	  MatrixXd		_alpha;
	  MatrixXd		 _T;	// system dipole field tensors
	  MatrixXd	 	_IDENT;	// a few temporaries for calculating eq 23
	  MatrixXd		_g;	
	  MatrixXd		_f;	
	  MatrixXd		_h;

	  MatR			_A;		// total system polarizability
	  VecR			_M;		// total system dipole moment (at time zero)

	  MatR_list		_vA;	// the tensor A for each timestep

	  bool	time_zero;

	  void CalculateTensors();
	  void CalculateTotalDipole ();
	  void CalculateLocalFieldCorrection ();
	  void CalculateTotalPolarizability ();

  };	// sfg-analyzer



  /*
  // Using the LAPACK solver for some simple systems
  extern "C" {

	void dgesv_(int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);

	void pdgesv_ (int *n, int *nrhs, double *A, int *ia, int *ja, int *desca, int* ipiv, double *B, int *ib, int *jb, int *descb, int *info); 

  } // extern
  */



}	// namespace morita

#endif
