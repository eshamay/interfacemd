#include "morita2002.h"

namespace morita {

  SFGAnalyzer::SFGAnalyzer (WaterSystemParams& wsp)
	: 
	  Analyzer<AmberSystem> (wsp),
	  _p(0), _alpha(0,0), _IDENT(0,0), _g(0,0), _f(0,0), time_zero(true)
  
  { return; }

  SFGAnalyzer::~SFGAnalyzer () {
  }

  void SFGAnalyzer::Setup () {

	if (time_zero) {
	  // load all the waters in the system
	  LoadWaters();

	  // all_wats will hold the collection of all waters in the system
	  for (Mol_it it = int_wats.begin(); it != int_wats.end(); it++) {
		MoritaH2O_ptr ptr (new MoritaH2O (*it));
		all_wats.push_back(ptr);
	  }
	}

	cutoff_wats.clear();
	std::copy(all_wats.begin(), all_wats.end(), back_inserter(cutoff_wats));

	// find the system center of mass
	VecR com (this->CenterOfMass<MoritaH2O_ptr>(all_wats));
	double cutoff = com[WaterSystem<AmberSystem>::axis];

	// only use waters found above a particular location in the slab (i.e. center of mass)
	std::pair<double,double> slice = std::make_pair(cutoff,WaterSystem<AmberSystem>::posmax);
	this->SliceWaters<MoritaH2O_ptr> (cutoff_wats, slice);

	int N = 3*cutoff_wats.size();

	_T.setZero(N,N);
	_p.setZero(N);
	_alpha.setZero(N,N);

  } // Setup


  void SFGAnalyzer::Analysis () {

	boost::timer t;
	// Sets up the p, alpha, and T tensors by calculating through each water's dipole moment, polarizability, and dipole field contributions
	this->CalculateTensors();

	// sums all the water dipole moments to get the total system dipole
	if (time_zero) {
	  this->CalculateTotalDipole();
	  time_zero = false;
	}

	t.restart();
	// determine the local field correction to each of the molecular polarizabilities
	this->CalculateLocalFieldCorrection ();
	printf ("local field correction:: ");
	std::cout << t.elapsed() << std::endl;

	/*
	// sum all the molecular polarizabilities to get the total system value
	this->CalculateTotalPolarizability ();
	*/

  } // Analysis


  void SFGAnalyzer::DataOutput () {

	/*

	// for now, only output the SSP and SPS components of the correlation function
	double a_sps, m_sps;
	double a_ssp, m_ssp;

	std::vector<double> sps, ssp;
	*/

	// *********** time-domain work *********** //
	/*
	for (tensor::tensor_it it = _vA.begin(); it != _vA.end(); it++) {
	  a_sps = ((*it)(0,1) + (*it)(2,1))/2.0;
	  a_ssp = ((*it)(0,0) + (*it)(0,2) + (*it)(2,0) + (*it)(2,2))/4.0;

	  m_sps = (_M(0) + _M(2))/2.0;
	  m_ssp = _M(1);
	  sps.push_back(a_sps*m_sps);
	  ssp.push_back(a_ssp*m_ssp);
	}
	*/


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

  void SFGAnalyzer::PostAnalysis () {
	return; 
  }


  void SFGAnalyzer::CalculateTensors() {
	// Calculate the dipole moment of each water, and then constructs the 3Nx3 tensor 'p'.
	if (time_zero)
	  std::for_each (cutoff_wats.begin(), cutoff_wats.end(), SetDipoleMoment());

	// Set up the polarizability (alpha) matrix similar to the method for the dipole moment
	std::for_each (cutoff_wats.begin(), cutoff_wats.end(), SetPolarizability());

	_T.setZero();
	for (unsigned int i = 0; i < cutoff_wats.size(); i++) {

	  if (time_zero) {
		_p.block(3*i,0,3,1) = cutoff_wats[i]->Dipole();
	  }

	  _alpha.block(3*i,3*i,3,3) = cutoff_wats[i]->Polarizability();

	  for (unsigned int j = i+1; j < cutoff_wats.size(); j++) {

		// Calculate the tensor 'T' which is formed of 3x3 matrix elements
		MatR dft = DipoleFieldTensor(cutoff_wats[i], cutoff_wats[j]);
		_T.block(3*i,3*j,3,3) = dft;
		_T.block(3*j,3*i,3,3) = dft;
	  }
	}



  }	// Calculate Tensors



  // take care of the polarization calculations for the first timestep
  void SFGAnalyzer::CalculateTotalDipole () {
	_M.Zero();
	for (Morita_it it = cutoff_wats.begin(); it != cutoff_wats.end(); it++)
	  _M += (*it)->Dipole();

  }	// calculate total dipole


  void SFGAnalyzer::CalculateLocalFieldCorrection () {
	// following equation 23 - (1 + T*alpha) f = h.
	// first here set up 1 + T*alpha, then solve the equation for f with a lapack routine

	// do lots of initialization before the polarizability calculations
	int N = 3*cutoff_wats.size();

	_g.setZero(N,N);

	char trans = 'N';
	double scale = 1.0;
	boost::timer t;
	t.restart();
	dgemm (&trans, &trans, &N, &N, &N, &scale, &_T(0,0), &N, &_alpha(0,0), &N, &scale, &_g(0,0), &N);
	std::cout << "matrix mult. " << t.elapsed() << std::endl;

	_IDENT.setIdentity(N,N);
	t.restart();
	_g += _IDENT;	// this is now 1 + T*alpha, a.k.a inverse of g
	std::cout << "matrix sum " << t.elapsed() << std::endl;

	// for now, f is 'h', the 3Nx3 block identity tensor
	_f.setZero(N,3);
	tensor::BlockIdentity(_f,3);

	int nrhs = 3;
	int ipiv[N];
	for (int i = 0; i < N; i++) ipiv[i] = 0;
	int info = 0;

	// now solve for f in the equation g*f = h
	t.restart();
	dgesv (&N, &nrhs, &_g(0,0), &N, ipiv, &_f(0,0), &N, &info);
	std::cout << "system solve:  " << t.elapsed() << std::endl;

  }	// Local field correction


  /*
  // calculates A for the current timestep
  void SFGAnalyzer::CalculateTotalPolarizability () {
	// determine 'A', the summed system polarizability.
	tensor::tensor_t fs (3,3);	// slices
	tensor::tensor_t as (3,3);
	tensor::tensor_t fa (3,3);

	_A.clear();
	for (int i = 0; i < cutoff_wats.size(); i++) {
	  // grab the piece of the local field tensor
	  fs.assign(project (_f, slice (3*i,1,3), slice(0,1,3)));
	  // and the alpha tensor
	  as.assign(project (_alpha, slice (3*i,1,3), slice(3*i,1,3)));
	  fa.assign(prod(as,fs));

	  _A.plus_assign(fa);
	}

	// push the polarizability tensor into the collection
	_vA.push_back (_A);
  }


  */
  //
  // Calculates the dipole moment vector of a given water molecule according to the method of Morita/Hynes 2002
  //
  void MoritaH2O::SetDipoleMoment () {

	this->SetAtoms();	// first get the bonds set in the water

	// determine the bondlength displacements
	dR1 = this->OH1()->Magnitude() - BONDLENGTH_EQ;
	dR1 *= sfg_units::ANG2BOHR;
	dR2 = this->OH2()->Magnitude() - BONDLENGTH_EQ;
	dR2 *= sfg_units::ANG2BOHR;

	// angle displacement from the equilibrium
	dA = (acos(this->Angle())*180.0/M_PI) - ANGLE_EQ;

	// symmetry-adapted coordinates
	X1 = dR1 + dR2;
	X2 = dR1 - dR2;

	// determine the charge from the fitting parameters
	qO = CHARGE_O_EQ + C1*X1 + C2*dA + C3*X1*X1 + C4*X2*X2 + C5*dA*dA + C6*X1*dA;
	double dq = C7*X2 + C8*X1*X2 + C9*X2*dA; // originally delta q1 - delta q2
	qH1 = (dq-qO)/2.0;
	qH2 = (-dq-qO)/2.0;

	// the dipole vector is thus classically determined (summing position * charge)
	this->_dipole.Zero();
	this->_dipole += (_o->Position() * qO);
	this->_dipole += (_h1->Position() * qH1);
	this->_dipole += (_h2->Position() * qH2);

  } // Set dipole moment




  void MoritaH2O::SetPolarizability () {
	_alpha1.Zero();

	_alpha1(0,0) = D1+D4*dR1;
	_alpha1(1,1) = D2+D5*dR1;
	_alpha1(2,2) = D3+D6*dR1+D7*dR1*dR1;

	_alpha2(0,0) = D1+D4*dR2;
	_alpha2(1,1) = D2+D5*dR2;
	_alpha2(2,2) = D3+D6*dR2+D7*dR2*dR2;

	_alpha.setZero();

	// These rotations will produce a polarizability tensor that is formed 
	// in the space-fixed frame (instead of the local or molecular frames
	// as written in the paper).
	this->DCMToLabMorita(z,1);
	_alpha = _alpha + (this->_DCM.transpose() * _alpha1 * this->_DCM);

	this->DCMToLabMorita(z,2);
	_alpha = _alpha + (this->_DCM.transpose() * _alpha2 * this->_DCM);
  }
  


  DipoleFieldTensor::DipoleFieldTensor (const MoritaH2O_ptr wat1, const MoritaH2O_ptr wat2)
  {
	VecR r = MDSystem::Distance(wat1->GetAtom("O"), wat2->GetAtom("O"));
	// work in atomic units (au)
	double distance = r.Magnitude() * sfg_units::ANG2BOHR;
	double ir3 = 1.0/pow(distance,3.0);
	double ir5 = 3.0/pow(distance,5.0);

	// calculate T as in eq. 10 of the Morita/Hynes 2002 paper
	for (unsigned i = 0; i < 3; i++) {
	  for (unsigned j = 0; j < 3; j++) {
		this->operator()(i,j) = ir3 - ir5*r[i]*r[j];
	  }
	}
  } //dipole field tensor c-tor


} // namespace morita



// used to calculate the SFG spectrum based on the morita/hynes 2002 method
int main (int argc, char **argv) {

  libconfig::Config cfg;
  cfg.readFile("system.cfg");

  std::string filename = cfg.lookup("analysis.morita2002.filename");
  libconfig::Setting &analysis = cfg.lookup("analysis");
  analysis.add("filename", libconfig::Setting::TypeString) = filename;

  WaterSystemParams wsp (cfg);

  morita::SFGAnalyzer sfg (wsp);

  sfg.SystemAnalysis ();

  return 0;
}
