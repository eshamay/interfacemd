#ifndef MORITA2002_H_
#define MORITA2002_H_

#include "../utility.h"
#include "../analysis.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <cmath>
#include <iostream>



//#include "../moritasfg2002.h"

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

	  MoritaH2O (const Molecule& molecule) : Water(molecule) { return; }
	  MoritaH2O (const Molecule * molecule) : Water(*molecule) { return; }

	  void SetDipoleMoment ();
	  void SetPolarizability ();

	protected:

	  VecR _dipole;
	  MatR _alpha1, _alpha2;	// polarizabilities of the two OH bonds
	  MatR _alpha;
	  MatR _DCM;

	  double dR1, dR2, dA;		// displacements of the OH bonds and the HOH angle from their equilibrium values
	  double X1, X2;
	  double qO, qH1, qH2;		// charges calculated from the MoritaHynes 2002 method
  };	// morita-h2o

  typedef std::vector<MoritaH2O *> Morita_ptr_vec;
  typedef Morita_ptr_vec::const_iterator Morita_it;



  namespace math {

	typedef matrix<double> matrix_t;
	typedef identity_matrix<double> id_matrix_t;

	template <typename T>
	class Tensor : public T {

	  public:
		Tensor (const int M, const int N) : T(M, N) { }
		Tensor (const matrix_t& t) : T(t) { }
		virtual ~Tensor() { }

		virtual void Print() const;

	}; // tensor

	typedef Tensor<matrix_t>	tensor_t;



	// A matrix that forms an Mx3 tensor where each 3x3 submatrix is the identity matrix
	class ColumnIdentityMatrix : public tensor_t {

	  public:
		ColumnIdentityMatrix (const int M) : tensor_t(M,3)
	  {
		id_matrix_t id (3,3);
		for (unsigned i = 0; i < size1()/3; i++)
		  project(*this, slice(3*i,1,3), slice(0,1,3)) = id;
	  }
	}; // Column Identity Matrix



	// dipole field tensor 'T' as used in the morita&hynes paper
	// it's a square 3Nx3N matrix, where N = number of particles
	class DipoleFieldTensor : public tensor_t {
	  public:
		DipoleFieldTensor (const Molecule* wat1, const Molecule* wat2);
	}; // Dipole field tensor

	typedef Tensor<symmetric_matrix<double, upper> > sym_tensor_t;



	class SymmetricMatrix : public sym_tensor_t {
	  public:
		SymmetricMatrix () : sym_tensor_t(0,0) { }
	}; // symmetric matrix


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
	  void PostAnalysis () { return; }

	private:
	  Morita_ptr_vec	_wats;
	  math::tensor_t	_p;
	  math::tensor_t	_alpha;
	  math::SymmetricMatrix _T;	// system dipole field tensors
  };	// sfg-analyzer



  namespace utilities {

	template <typename T, typename U>
	  class ConvertPointerType : public std::unary_function<T,U> {
		public:
		  U operator() (T t) const { return static_cast<U>(t); }
	  };

	template <typename Iter1, typename Iter2>
	  void ConvertContainerElementTypes (Iter1 pBegin, Iter1 pEnd, Iter2 pBegin2)
	  {
		typedef typename std::iterator_traits<Iter1>::value_type value_t1;
		typedef typename std::iterator_traits<Iter2>::value_type value_t2;

		std::transform(pBegin, pEnd, pBegin2, ConvertPointerType<value_t1, value_t2>());
	  }

	template <typename T, typename U>
	  class ConvertAClass : public std::unary_function<T *,U *> {
		public:
		  U * operator() (T * t) const {
			return new U(t);
		  }
	  };


  } // namespace utilities



}	// namespace morita

#endif
