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

	template <class U>
	class SFGAnalyzer : public Analyzer<U> {
	  public:
		SFGAnalyzer (WaterSystemParams& params);
		~SFGAnalyzer ();

		virtual void Setup ();
		virtual void Analysis ();
		virtual void DataOutput ();
		void PostAnalysis () { return; }

	  protected:
		Morita_ptr_vec		all_wats;
		Morita_ptr_vec		analysis_wats;
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

		virtual void SelectAnalysisWaters () = 0;	// routine that determines which waters will be used in the analysis

		void CalculateTensors();
		virtual void SetAnalysisWaterDipoleMoments () = 0;

		void CalculateTotalDipole ();
		void CalculateTotalPolarizability ();
		void CalculateLocalFieldCorrection ();

		sfgdata_pair_t SSP_SPS_Result (const int s1, const int s2, const int p, const MatR& a) const;
	};	// sfg-analyzer


  template <class U>
	SFGAnalyzer<U>::SFGAnalyzer (WaterSystemParams& wsp)
	: 
	  Analyzer<U> (wsp),
	  _p(0), _alpha(0,0), _IDENT(0,0), _g(0,0), _f(0,0), time_zero(true)

  { return; }

  template <class U>
	SFGAnalyzer<U>::~SFGAnalyzer () {
	}

  template <class U>
	void SFGAnalyzer<U>::Setup () {

	  if (time_zero) {
		// load all the waters in the system
		this->LoadWaters();

		// all_wats will hold the collection of all waters in the system
		for (Mol_it it = this->int_wats.begin(); it != this->int_wats.end(); it++) {
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
	  _p.setZero(N);
	  _alpha.setZero(N,N);

	} // Setup



  template <class U>
	void SFGAnalyzer<U>::Analysis () {

	  // Sets up the p, alpha, and T tensors by calculating through each water's dipole moment, polarizability, and dipole field contributions
	  this->CalculateTensors();

	  // sums all the water dipole moments to get the total system dipole
	  //if (time_zero) {
		this->CalculateTotalDipole();
		time_zero = false;
	 // }

	  //t.restart();
	  // determine the local field correction to each of the molecular polarizabilities
	  this->CalculateLocalFieldCorrection ();
	  //std::cout << "local field correction:: " << t.elapsed() << std::endl;

	  // sum all the molecular polarizabilities to get the total system value
	  this->CalculateTotalPolarizability ();

	} // Analysis


  // calculates the ssp and sps value of the autocorrelation function for a given set of pqr space-frame coordinates
  template <class U>
	sfgdata_pair_t SFGAnalyzer<U>::SSP_SPS_Result (const int s1, const int s2, const int p, const MatR& a) const {
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


  template <class U>
	void SFGAnalyzer<U>::DataOutput () {

	  rewind (this->output);

	  // try this - for each timestep, output the vector dipole, and the matrix polarizability of the system
	  for (MatR_it it = _vA.begin(); it != _vA.end(); it++) {

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

	  fflush(this->output);


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

	} // Data Output

  template <class U>
	void SFGAnalyzer<U>::CalculateTensors() {
	  // Calculate the dipole moment of each water, and then constructs the 3Nx3 tensor 'p'.
	  if (time_zero)
		this->SetAnalysisWaterDipoleMoments ();

	  // Set up the polarizability (alpha) matrix similar to the method for the dipole moment
	  std::for_each (analysis_wats.begin(), analysis_wats.end(), std::mem_fun(&MoritaH2O::SetPolarizability));

	  _T.setZero();
	  for (unsigned int i = 0; i < analysis_wats.size(); i++) {

		if (time_zero) {
		  _p.block(3*i,0,3,1) = analysis_wats[i]->Dipole();
		}

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
	void SFGAnalyzer<U>::CalculateTotalDipole () {

	  VecR_vec dipoles;
	  std::transform (analysis_wats.begin(), analysis_wats.end(), std::back_inserter(dipoles), std::mem_fun<VecR,Molecule>(&Molecule::Dipole));

	  _M = std::accumulate (dipoles.begin(), dipoles.end(), VecR());

	}	// calculate total dipole


  template <class U>
	void SFGAnalyzer<U>::CalculateLocalFieldCorrection () {
	  // following equation 23 - (1 + T*alpha) f = h.
	  // first here set up 1 + T*alpha, then solve the equation for f with a lapack routine

	  // do lots of initialization before the polarizability calculations
	  int N = 3*analysis_wats.size();

	  _g.setZero(N,N);

	  char trans = 'N';
	  double scale = 1.0;
	  //boost::timer t;
	  //t.restart();
	  dgemm (&trans, &trans, &N, &N, &N, &scale, &_T(0,0), &N, &_alpha(0,0), &N, &scale, &_g(0,0), &N);
	  //std::cout << "blas matrix mult. " << t.elapsed() << std::endl;

	  _IDENT.setIdentity(N,N);

	  //t.restart();
	  _g += _IDENT;	// this is now 1 + T*alpha, a.k.a inverse of g
	  //std::cout << "eigen matrix sum " << t.elapsed() << std::endl;

	  // for now, f is 'h', the 3Nx3 block identity tensor
	  _h.setZero(N,3);
	  tensor::BlockIdentity(_h,3);

	  // now solve for f in the equation g*f = h 
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
	void SFGAnalyzer<U>::CalculateTotalPolarizability () {
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
