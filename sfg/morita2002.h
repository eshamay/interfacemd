
#define DEBUG__	1
#ifndef MORITA2002_H_
#define MORITA2002_H_

#define MPI__	1

//#include "../utility.h"
#include "../analysis.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <iostream>


namespace morita {

  using namespace boost::numeric::ublas;
  /*	From the morita 2002 water model
		const double BONDLENGTH_EQ	= 0.9575;		// in angstroms
		const double ANGLE_EQ 		= 104.51*M_PI/180.0;		// in radians
		const double CHARGE_H_EQ		= 0.3285;		// charge units (atomic units?)
		const double CHARGE_O_EQ		= -0.6570;
   */

  /* from the equilibrium values of the spcfw model */
  const double BONDLENGTH_EQ	= 1.012;		// in angstroms
  const double ANGLE_EQ 		= 113.24*M_PI/180.0;		// in radians
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

  typedef std::vector<MoritaH2O *> Morita_ptr_vec;
  typedef Morita_ptr_vec::const_iterator Morita_it;


  namespace math {


	// dipole field tensor 'T' as used in the morita&hynes paper
	// it's a square 3Nx3N matrix, where N = number of particles
	class DipoleFieldTensor : public tensor::tensor_t {
	  public:
		DipoleFieldTensor (const Molecule* wat1, const Molecule* wat2);
	}; // Dipole field tensor


  } // math



  class SetDipoleMoment : public std::unary_function<MoritaH2O *,void> {
	public:
	  void operator() (MoritaH2O * wat) { wat->SetDipoleMoment(); }
  };

  class SetPolarizability : public std::unary_function<MoritaH2O *,void> {
	public:
	  void operator() (MoritaH2O * wat) { wat->SetPolarizability(); }
  };



  class SFGAnalyzer : public Analyzer<AmberSystem> {
	public:
	  SFGAnalyzer (WaterSystemParams& params);
	  ~SFGAnalyzer ();

	  void Setup ();
	  void Analysis ();
	  void DataOutput (const unsigned int timestep);
	  void PostAnalysis ();

	private:

	  Morita_ptr_vec		_wats;
	  vector_t				_p;
	  tensor::tensor_t		_alpha;
	  tensor::SymmetricMatrix _T;	// system dipole field tensors
	  tensor::id_matrix_t 	_IDENT;	// a few temporaries for calculating eq 23
	  tensor::tensor_t		_Talpha;
	  tensor::tensor_t		_ginv;	
	  tensor::tensor_t		_h;	
	  tensor::tensor_t		_f;	

  };	// sfg-analyzer



  namespace utilities {

	template <typename T>
	  class DeletePointer : public std::unary_function<T,void> {
		public:
		  void operator() (T t) { delete t; }
	  };  // delete pointer


	template <typename T, typename U>
	  class ConvertPointerType : public std::unary_function<T,U> {
		public:
		  U operator() (T t) const { return static_cast<U>(t); }
	  };  // convert pointer type


	template <typename Iter1, typename Iter2>
	  void ConvertContainerElementTypes (Iter1 pBegin, Iter1 pEnd, Iter2 pBegin2)
	  {
		typedef typename std::iterator_traits<Iter1>::value_type value_t1;
		typedef typename std::iterator_traits<Iter2>::value_type value_t2;

		std::transform(pBegin, pEnd, pBegin2, ConvertPointerType<value_t1, value_t2>());
	  } // convert container element types


	template <typename T, typename U>
	  class MakeDerivedFromPointer : public std::unary_function<T *,U *> {
		public:
		  U * operator() (T * t) const {
			return new U(*t);
		  }
	  };  // make derived from pointer


  } // namespace utilities


  // Using the LAPACK solver for some simple systems
  extern "C" {

	void dgesv_(int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);

	void pdgesv_ (int *n, int *nrhs, double *A, int *ia, int *ja, int *desca, int* ipiv, double *B, int *ib, int *jb, int *descb, int *info); 

  } // extern



}	// namespace morita

#endif
