#ifndef MORITA2002_H_
#define MORITA2002_H_

#include "../analysis.h"
#include "sfgunits.h"
#include "moritah2o.h"
#include "../tensor.h"
#include <Eigen/LU>
#include <iostream>
//#include <fftw3.h>
#include <mkl_blas.h>
#include <mkl_lapack.h>


namespace morita {

  typedef std::pair<double,double> sfgdata_pair_t;

  USING_PART_OF_NAMESPACE_EIGEN

	/*!
	 * A pure virtual base class for performing an SFG analysis based on the method of Morita/Hynes (J. Phys. Chem. B 2002, 106, 673-685). Depending on the specific MD system used (Amber, CP2K, etc) the pure virtual methods define the actions to be taken to alter how various components are defined, and customize the analysis.
	 */
  template <class U>
  class Morita2002Analysis : public AnalysisSet< Analyzer<U> > {
	public:
	  Morita2002Analysis (std::string, std::string);
	  ~Morita2002Analysis ();

	  typedef Analyzer<U> system_t;

	  virtual void Setup (system_t&);
	  virtual void Analysis (system_t&);
	  virtual void DataOutput (system_t&);

	protected:
	  //! all the waters in the MD system
	  Morita_ptr_vec		all_wats;		
	  //! Water molecules that will be analyzed by the morita 2002 routines
	  Morita_ptr_vec		analysis_wats;	
	  //VectorXd				_p;
	  //! The sum-total of all molecular polarizabilities of the water molecules in the system
	  MatrixXd		_alpha;
	  //! System dipole field tensor
	  MatrixXd		 _T;	
	  //! Identity matrix used in calculation of eq.23 of the Morita/Hynes 2002 method
	  MatrixXd	 	_IDENT;	
	  MatrixXd		_g;	
	  //! Local field correction tensor
	  MatrixXd		_f;	
	  //! Block-identity matrix
	  MatrixXd		_h;
	  
	  //! total system polarizability
	  MatR			_A;		
	  //! total system dipole moment (at time zero)
	  VecR			_M;		

	  //! Collection of all total system dipole moments for each timestep
	  VecR_list		_vM;	
	  //! the total system polarizability matrix for each timestep
	  MatR_list		_vA;	

	  bool	time_zero;

	  /*!
	   * Determines which of the molecules (waters) in the system are to be used for the ensuing analysis. Only waters that have been copied into the analysis_wats container will be used.
	   */
	  virtual void SelectAnalysisWaters () = 0;	// routine that determines which waters will be used in the analysis

	  /*!
	   * Method of calculating the various tensors (molecular dipole moments, polarizabilities, dipole field tensor, etc) used in the ensuing analysis
	   */
	  void CalculateTensors();

	  //! Sets the dipole moment of each water used in the analysis.
	  virtual void SetAnalysisWaterDipoleMoments () = 0;

	  void CalculateTotalDipole ();
	  void CalculateTotalPolarizability ();
	  void CalculateLocalFieldCorrection ();

	  // sfgdata_pair_t SSP_SPS_Result (const int s1, const int s2, const int p, const MatR& a) const;
  };	// sfg-analyzer


  template <class U>
	Morita2002Analysis<U>::Morita2002Analysis (std::string desc, std::string fn)
	: 
	  AnalysisSet< Analyzer<U> > (desc, fn),
	  //_p(0), 
	  _alpha(0,0), _IDENT(0,0), _g(0,0), _f(0,0), time_zero(true) 
  { return; }



  template <class U> Morita2002Analysis<U>::~Morita2002Analysis () { } 


  template <class U>
	void Morita2002Analysis<U>::Setup (system_t& t) {


	  if (time_zero) {
		AnalysisSet < Analyzer<U> >::Setup(t);

		// load all the waters in the system
		t.LoadWaters();

		// all_wats will hold the collection of all waters in the system
		for (Mol_it it = t.int_wats.begin(); it != t.int_wats.end(); it++) {
		  MoritaH2O_ptr ptr (new MoritaH2O (*it));
		  all_wats.push_back(ptr);
		}

	  }

	  // analysis_wats will hold only those waters that will by analyzed
	  analysis_wats.clear();
	  std::copy(all_wats.begin(), all_wats.end(), back_inserter(analysis_wats));

	  // here filter out the waters to use for analysis
	  this->SelectAnalysisWaters ();

	  int N = 3*analysis_wats.size();
	  _T.setZero(N,N);
	  //_p.setZero(N);
	  _alpha.setZero(N,N);

	} // Setup



  template <class U>
	void Morita2002Analysis<U>::Analysis (system_t& t) {

	  // Sets up the p, alpha, and T tensors by calculating through each water's dipole moment, polarizability, and dipole field contributions
	  this->CalculateTensors();

	  // determine the local field correction to each of the molecular polarizabilities
	  this->CalculateLocalFieldCorrection ();

	  // sum all the molecular polarizabilities to get the total system value
	  this->CalculateTotalPolarizability ();



	  for (unsigned int i = 0; i < 3; i++) {
		fprintf (t.Output(), "%8.3f ", _M(i));
	  } 
	  // output in row-major order
	  for (unsigned int i = 0; i < 3; i++) {
		for (unsigned int j = 0; j < 3; j++) {
		  fprintf (t.Output(), "% 8.3f", _A(i,j));
		}
	  }
	  fprintf (t.Output(),"\n");

	} // Analysis


  // calculates the ssp and sps value of the autocorrelation function for a given set of pqr space-frame coordinates
  /*
  template <class U>
	sfgdata_pair_t Morita2002Analysis<U>::SSP_SPS_Result (const int s1, const int s2, const int p, const MatR& a) const {
	  // for now, only output the SSP and SPS components of the correlation function
	  double a_sps, m_sps;
	  double a_ssp, m_ssp;

	  a_sps = a(s1,p);
	  m_sps = _M(s1);

	  a_ssp = a(s1,s1);
	  m_ssp = _M(p);
	  //a_sps = (a(s1,p) + a(s2,p))/2.0;
	  //a_ssp = (a(s1,s1) + a(s2,s2))/2.0;

	  //m_sps = (_M(s1) + _M(s2))/2.0;
	  //m_ssp = _M(p);

	  return std::make_pair (a_ssp*m_ssp, a_sps*m_sps);
	}
	*/


  template <class U>
  void Morita2002Analysis<U>::DataOutput (system_t& t) {


  // try this - for each timestep, output the vector dipole, and the matrix polarizability of the system

  // *********** time-domain output work *********** //
  // Output will be 6 columns. 2 each for the X,Y, and Z axes as the "P" axis, and each pair of columns will have ssp, sps data
  /*
	 sfgdata_pair_t datapair;
	 for (MatR_it it = _vA.begin(); it != _vA.end(); it++) {
  //datapair = SSP_SPS_Result(1,2,0,*it);	// for the X-axis
  //fprintf (this->output, "% 12.6f % 12.6f", datapair.first, datapair.second);
  //datapair = SSP_SPS_Result(0,2,1,*it);	// for the Y-axis
  //fprintf (this->output, "% 12.6f % 12.6f", datapair.first, datapair.second);
  datapair = SSP_SPS_Result(1,2,0,*it);	// for the Z-axis
  fprintf (this->output, "% 12.6f % 12.6f", datapair.first, datapair.second);
  datapair = SSP_SPS_Result(2,1,0,*it);	// for the Z-axis
  fprintf (this->output, "% 12.6f % 12.6f", datapair.first, datapair.second);
  datapair = SSP_SPS_Result(0,2,1,*it);	// for the Z-axis
  fprintf (this->output, "% 12.6f % 12.6f", datapair.first, datapair.second);
  datapair = SSP_SPS_Result(2,0,1,*it);	// for the Z-axis
  fprintf (this->output, "% 12.6f % 12.6f", datapair.first, datapair.second);
  datapair = SSP_SPS_Result(0,1,2,*it);	// for the Z-axis
  fprintf (this->output, "% 12.6f % 12.6f", datapair.first, datapair.second);
  datapair = SSP_SPS_Result(1,0,2,*it);	// for the Z-axis
  fprintf (this->output, "% 12.6f % 12.6f\n", datapair.first, datapair.second);
  }
  */

	fflush(t.Output());
  }


  // ********  fftw work ******** //
  /*
	 double *in;
	 fftw_complex *sps_fft, *ssp_fft;
	 fftw_plan p;

	 in = (double*) fftw_malloc(sizeof(double)*sps.size());

  // sps transform
  std::copy(sps.begin(), sps.end(), in);
  sps_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*sps.size());
  p = fftw_plan_dft_r2c_1d(sps.size(),in,sps_fft,FFTW_DESTROY_INPUT);
  fftw_execute(p);

  // ssp transform
  std::copy(ssp.begin(), ssp.end(), in);
  ssp_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ssp.size());
  p = fftw_plan_dft_r2c_1d(ssp.size(),in,ssp_fft,FFTW_DESTROY_INPUT);
  fftw_execute(p);

  // output one column each for
  // time-domain data, real, imaginary, magnitude
  rewind(output);
  for (int i = 0; i < sps.size(); i++) {
  fprintf (output, "% 23.8f % 23.8f % 23.8f % 23.8f % 23.8f % 23.8f % 23.8f % 23.8f\n",
  sps[i], sps_fft[i][0], sps_fft[i][1], sqrt(sps_fft[i][0]*sps_fft[i][0] + sps_fft[i][1]*sps_fft[i][1]),
  ssp[i], ssp_fft[i][0], ssp_fft[i][1], sqrt(ssp_fft[i][0]*ssp_fft[i][0] + ssp_fft[i][1]*ssp_fft[i][1])
  );
  }

  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(sps_fft);
  fftw_free(ssp_fft);


  if (!(timestep == timesteps))
  this->Setup();
  */

  //} // Data Output

  template <class U>
	void Morita2002Analysis<U>::CalculateTensors() {

	  // Calculate the dipole moment of each water, and then constructs the 3Nx3 tensor 'p'.
	  //if (time_zero)
	  this->SetAnalysisWaterDipoleMoments ();
	  // sums all the water dipole moments to get the total system dipole

	  this->CalculateTotalDipole();

	  // Set up the polarizability (alpha) matrix similar to the method for the dipole moment
	  std::for_each (analysis_wats.begin(), analysis_wats.end(), std::mem_fun(&MoritaH2O::SetPolarizability));

	  _T.setZero();
	  for (unsigned int i = 0; i < analysis_wats.size(); i++) {


		//if (time_zero) {
		//_p.block(3*i,0,3,1) = analysis_wats[i]->Dipole();
		//}

		_alpha.block(3*i,3*i,3,3) = analysis_wats[i]->Polarizability();

		for (unsigned int j = i+1; j < analysis_wats.size(); j++) {

		  // Calculate the tensor 'T' which is formed of 3x3 matrix elements
		  MatR dft = DipoleFieldTensor(analysis_wats[i], analysis_wats[j]);
		  _T.block(3*i,3*j,3,3) = dft;
		  _T.block(3*j,3*i,3,3) = dft;
		}
	  }

	}	// Calculate Tensors



  // take care of the polarization calculations for the first timestep
  template <class U>
	void Morita2002Analysis<U>::CalculateTotalDipole () {
	  VecR_vec dipoles;
	  std::transform (analysis_wats.begin(), analysis_wats.end(), std::back_inserter(dipoles), std::mem_fun<VecR,Molecule>(&Molecule::Dipole));

	  // sum all the molecular dipoles together to get the total system dipole
	  _M = std::accumulate (dipoles.begin(), dipoles.end(), VecR());

	  // in case it's needed - keep the dipole for every timestep in a running list
	  _vM.push_back(_M);

	}	// calculate total dipole


  template <class U>
	void Morita2002Analysis<U>::CalculateLocalFieldCorrection () {
	  // following equation 23 - (1 + T*alpha) f = h.
	  // first, set up _g = (1 + T*alpha)
	  // then solve the equation for f with a lapack routine
	  //
	  // note: g = inv(1+T*alpha)
	  // so f = g*h --- solved with dsgesv from lapack

	  // do lots of initialization before the polarizability calculations
	  int N = 3*analysis_wats.size();

	  _g.setZero(N,N);

	  char trans = 'N';
	  double scale = 1.0;
	  //boost::timer t;
	  //t.restart();
	  // set g = T*alpha
	  dgemm (&trans, &trans, &N, &N, &N, &scale, &_T(0,0), &N, &_alpha(0,0), &N, &scale, &_g(0,0), &N);
	  //std::cout << "blas matrix mult. " << t.elapsed() << std::endl;

	  _IDENT.setIdentity(N,N);

	  //t.restart();
	  // now g = 1+T*alpha
	  _g += _IDENT;
	  //std::cout << "eigen matrix sum " << t.elapsed() << std::endl;

	  // for now, f is 'h', the 3Nx3 block identity tensor
	  _h.setZero(N,3);
	  tensor::BlockIdentity(_h,3);

	  // now solve for f in the equation g*f = h using the lapack dsgesv
	  int nrhs = 3;
	  int ipiv[N];
	  for (int i = 0; i < N; i++) ipiv[i] = 0;
	  int info = 0;

	  double work[N*nrhs];
	  float swork[N*(N+nrhs)];
	  int iter;
	  _f.setZero(N,3);

	  //t.restart();
	  dsgesv (&N, &nrhs, &_g(0,0), &N, ipiv, &_h(0,0), &N, &_f(0,0), &N, work, swork, &iter, &info);
	  if (info != 0){
		std::cout << "DSGESV.info parameter had a value of " << info << " meaning that something went wrong!" << std::endl;
		exit(1);
	  }
	  //std::cout << "DSGESV (iterative) system solve:  " << t.elapsed() << std::endl;

	  //t.restart();
	  //dgesv (&N, &nrhs, &_g(0,0), &N, ipiv, &_f(0,0), &N, &info);
	  //std::cout << "DGESV system solve:  " << t.elapsed() << std::endl;

	  //MatrixXd _x;
	  //t.restart();
	  //_g.lu().solve(_f, &_x);
	  //std::cout << "Eigen LU solve:  " << t.elapsed() << std::endl;

	}	// Local field correction


  // calculates A for the current timestep
  template <class U>
	void Morita2002Analysis<U>::CalculateTotalPolarizability () {
	  // determine 'A', the summed (total) system polarizability.
	  Matrix3d fa (3,3);

	  int N = analysis_wats.size();
	  _A.setZero(N,N);
	  for (int i = 0; i < N; i++) {
		// grab the piece of the local field tensor and find the inner product with alpha
		fa = _alpha.block(3*i,3*i,3,3) * _f.block(3*i,0,3,3);
		_A += fa;
	  }
	  // push the polarizability tensor into the collection
	  _vA.push_back (_A);
	}	// Calculate total polarizability


  /*
  // Using the LAPACK solver for some simple systems
  extern "C" {

  void dgesv_(int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);

  void pdgesv_ (int *n, int *nrhs, double *A, int *ia, int *ja, int *desca, int* ipiv, double *B, int *ib, int *jb, int *descb, int *info); 

  } // extern
  */


}	// namespace morita

#endif
