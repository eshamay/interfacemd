#include "morita2002.h"

namespace morita {

  using namespace boost::numeric::ublas;

  SFGAnalyzer::SFGAnalyzer (WaterSystemParams& wsp)
	: 
	  Analyzer<AmberSystem> (wsp),
	  _p(0), _alpha(0,0), _IDENT(0), _g(0,0), _f(0,0), _A(3,3), time_zero(true)
  
  { return; }

  SFGAnalyzer::~SFGAnalyzer () {
	for (Morita_it it = _wats.begin(); it != _wats.end(); it++) { delete (*it); }
  }

  void SFGAnalyzer::Setup () {

	// load all the waters in the system
	LoadWaters();

	// clear out and then reload the new set of waters in the system
	std::for_each(_wats.begin(), _wats.end(), utilities::DeletePointer<MoritaH2O*>());
	_wats.clear();
	std::transform(int_wats.begin(), int_wats.end(), std::back_inserter(_wats), utilities::MakeDerivedFromPointer<Molecule, MoritaH2O>());

	int N = _wats.size();

	_T.resize(3*N); _T.clear();
	_p.resize(3*N); _p.clear();
	_alpha.resize(3*N, 3*N); _alpha.clear();

  } // Setup

  void SFGAnalyzer::Analysis () {

	_T.clear();

	// Sets up the p, alpha, and T tensors by calculating through each water's dipole moment, polarizability, and dipole field contributions
	this->CalculateTensors();

	// sums all the water dipole moments to get the total system dipole
	if (time_zero)
	  this->CalculateTotalDipole();

	// determine the local field correction to each of the molecular polarizabilities
	this->CalculateLocalFieldCorrection ();

	// sum all the molecular polarizabilities to get the total system value
	this->CalculateTotalPolarizability ();

  } // Analysis


  void SFGAnalyzer::DataOutput (const unsigned int timestep) {

  } // Data Output




  void SFGAnalyzer::CalculateTensors() {
	// Calculate the dipole moment of each water, and then constructs the 3Nx3 tensor 'p'.
	if (time_zero)
	  std::for_each (_wats.begin(), _wats.end(), SetDipoleMoment());

	// Set up the polarizability (alpha) matrix similar to the method for the dipole moment
	std::for_each (_wats.begin(), _wats.end(), SetPolarizability());

	for (unsigned int i = 0; i < _wats.size(); i++) {

	  if (time_zero)
		project(_p, slice(3*i,1,3)) = _wats[i]->Dipole();

	  project(_alpha, slice(3*i,1,3), slice(3*i,1,3)) = _wats[i]->Polarizability();

	  for (unsigned int j = i+1; j < _wats.size(); j++) {

		if (i == j) continue;

		// Calculate the tensor 'T' which is formed of 3x3 matrix elements
		project(_T, slice(3*i,1,3), slice(3*j,1,3)) = math::DipoleFieldTensor(_wats[i], _wats[j]);
	  }
	}
  }



  // take care of the polarization calculations for the first timestep
  void SFGAnalyzer::CalculateTotalDipole () {
	_M.Zero();
	for (Morita_it it = _wats.begin(); it != _wats.end(); it++)
	  _M.plus_assign((*it)->Dipole());

	time_zero = false;

  }	// calculate total dipole


  void SFGAnalyzer::CalculateLocalFieldCorrection () {
	// following equation 23 - (1 + T*alpha) f = h.
	// first here set up 1 + T*alpha, then solve the equation for f with a lapack routine

	// do lots of initialization before the polarizability calculations
	int N = 3*_wats.size();

	_IDENT.resize(N);
	_g.resize(N,N);	_g.clear();

	char trans = 'N';
	double scale = 1.0;
	blas::dgemm_ (&trans, &trans, &N, &N, &N, &scale, &_T(0,0), &N, &_alpha(0,0), &N, &scale, &_g(0,0), &N);

	_g.plus_assign (_IDENT);	// this is now 1 + T*alpha, a.k.a inverse of g

	// for now, f is 'h', the 3Nx3 block identity tensor
	_f.resize(N,3); _f.clear();
	tensor::tensor_t::BlockIdentity(_f,3);

	int nrhs = 3;
	int ipiv[N];
	int info;

	// now solve for f in the equation g*f = h
	lapack::dgesv_ (&N, &nrhs, &_g(0,0), &N, ipiv, &_f(0,0), &N, &info);

  }	// Local field correction


  // calculates A for the current timestep
  void SFGAnalyzer::CalculateTotalPolarizability () {
	// determine 'A', the summed system polarizability.
	tensor::tensor_t fs (3,3);	// slices
	tensor::tensor_t as (3,3);
	tensor::tensor_t fa (3,3);

	_A.clear();
	for (int i = 0; i < _wats.size(); i++) {
	  // grab the piece of the local field tensor
	  fs.assign(project (_f, slice (3*i,1,3), slice(0,1,3)));
	  // and the alpha tensor
	  as.assign(project (_alpha, slice (3*i,1,3), slice(3*i,1,3)));
	  fa.assign(prod(as,fs));

	  _A.plus_assign(fa);
	}
  }




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
	this->_dipole += _o->Position() * qO;
	this->_dipole += _h1->Position() * qH1;
	this->_dipole += _h2->Position() * qH2;
  } // Set dipole moment




  void MoritaH2O::SetPolarizability () {
	_alpha1.Zero();

	_alpha1.Set(0,0, D1+D4*dR1);
	_alpha1.Set(1,1, D2+D5*dR1);
	_alpha1.Set(2,2, D3+D6*dR1+D7*dR1*dR1);

	_alpha2.Set(0,0, D1+D4*dR2);
	_alpha2.Set(1,1, D2+D5*dR2);
	_alpha2.Set(2,2, D3+D6*dR2+D7*dR2*dR2);

	_alpha.Zero();

	// These rotations will produce a polarizability tensor that is formed 
	// in the space-fixed frame (instead of the local or molecular frames
	// as written in the paper).
	this->DCMToLabMorita(z,1);
	_alpha = _alpha + (this->_DCM.Transpose() * _alpha1 * this->_DCM);

	this->DCMToLabMorita(z,2);
	_alpha = _alpha + (this->_DCM.Transpose() * _alpha2 * this->_DCM);
  }



  namespace math {


	DipoleFieldTensor::DipoleFieldTensor (const Molecule* wat1, const Molecule* wat2)
	  :
		tensor::tensor_t(3,3)
	{
	  VecR r = MDSystem::Distance (wat1->GetAtom("O"), wat2->GetAtom("O"));
	  // work in atomic units (au)
	  double distance = r.Magnitude() * sfg_units::ANG2BOHR;
	  double ir3 = 1.0/pow(distance,3.0);
	  double ir5 = 3.0/pow(distance,5.0);

	  // calculate T as in eq. 10 of the Morita/Hynes 2002 paper
	  for (unsigned i = 0; i < 3; i++) {
		for (unsigned j = 0; j < 3; j++) {
		  (*this)(i,j) = ir3 - ir5*r[i]*r[j];
		}
	  }
	} //DFT::DFT


  }	// namespace math

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
