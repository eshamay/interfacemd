#include "morita2002.h"

namespace morita {

  using namespace boost::numeric::ublas;
  namespace bnu = boost::numeric::ublas;

  SFGAnalyzer::SFGAnalyzer (WaterSystemParams& wsp)
	: 
	  Analyzer<AmberSystem> (wsp),
	  _p(0), _alpha(0,0), _IDENT(0), _Talpha(0,0), _ginv(0,0), _h(0,0), _f(0,0)
  { return; }

  SFGAnalyzer::~SFGAnalyzer () {
  }




  void SFGAnalyzer::Setup () {

	int N;

	if (MPIStatus()) {
	  // load all the waters in the system
	  LoadWaters();

	  // clear out and then reload the new set of waters in the system
	  std::for_each(_wats.begin(), _wats.end(), utilities::DeletePointer<MoritaH2O*>());
	  _wats.resize(int_wats.size());
	  std::transform(int_wats.begin(), int_wats.end(), _wats.begin(), utilities::MakeDerivedFromPointer<Molecule, MoritaH2O>());

	  //int N = _wats.size();
	  //cout << "using " << N << " waters" << endl;
	  N = _wats.size();
	  _wats.resize(200);
	  N = _wats.size();
	}


	if (wsp.mpisys)
	  boost::mpi::broadcast(wsp.mpi_world,N,0);

	_T.resize(3*N); _T.clear();
	_p.resize(3*N); _p.clear();
	_alpha.resize(3*N, 3*N); _alpha.clear();


  } // Setup

  void SFGAnalyzer::Analysis () {

	int N;
	if (MPIStatus()) {
	  // Calculate the tensor 'T' which is formed of 3x3 matrix elements
	  _T.clear();
	  for (unsigned int i = 0; i < _wats.size()/3; i++) {
		for (unsigned int j = i+1; j < _wats.size()/3; j++) {
		  if (i == j) continue;
		  project(_T, bnu::slice(3*i,1,3), bnu::slice(3*j,1,3)) = math::DipoleFieldTensor(_wats[i], _wats[j]);
		}
	  }

	  // Calculate the dipole moment of each water, and then constructs the 3Nx3 tensor 'p'.
	  std::for_each (_wats.begin(), _wats.end(), SetDipoleMoment());
	  for (unsigned i = 0; i < _wats.size()/3; i++){
		project(_p, bnu::slice(3*i,1,3)) = _wats[i]->Dipole();
	  }

	  /*
		 for (int i = 0; i < 3; i++)
		 _wats[i]->Dipole().Print();

		 cout << endl;
		 vector_slice<vector_t> s (_p, slice(0,1,9));
		 for (int i = 0; i < 9; i++)
		 cout << s(i) << endl;
	   */


	  // Set up the polarizability (alpha) matrix similar to the method for the dipole moment
	  std::for_each (_wats.begin(), _wats.end(), SetPolarizability());
	  for (unsigned int i = 0; i < _wats.size()/3; i++) {
		project(_alpha, bnu::slice(3*i,1,3), bnu::slice(3*i,1,3)) = _wats[i]->Polarizability();
	  }

	  // following equation 23 - (1 + T*alpha) f = h.
	  // first here set up 1 + T*alpha

	  // do lots of initialization before the calculations
	  N = 3*_wats.size();
	}

	if (wsp.mpisys)
	  boost::mpi::broadcast(wsp.mpi_world,N,0);

	_IDENT.resize(N);
	_Talpha.resize(N,N); _Talpha.clear();
	_ginv.resize(N,N);	_ginv.clear();
	_f.resize(N,3);
	_h.resize(N,3);
	tensor::tensor_t::BlockIdentity(_h,3);

	if (MPIStatus())
	  _Talpha.assign (prod(_T, _alpha));
	_ginv.assign (_IDENT);
	_ginv.plus_assign (_Talpha);	// this is now 1 + T*alpha, a.k.a inverse of g

	// a VERY long operation! avoid this :)
	//tensor::tensor_t::Inverse(_ginv, _g);

	int_vector_t Pivot(N);
	int info;
	int NRHS = 3;
	//_f.assign(_h);
	//dgesv_ (&N,&NRHS,&_ginv(0,0),&N,&Pivot(0),&_h(0,0),&N,&info);
	//pdgesv_ (&N, &NRHS, );

	if (MPIStatus())
	  _h.Print();


  } // Analysis



  void SFGAnalyzer::DataOutput (const unsigned int timestep) {

  }

  void SFGAnalyzer::PostAnalysis () {
	if (MPIStatus()) {
	  std::for_each(_wats.begin(), _wats.end(), utilities::DeletePointer<MoritaH2O*>());
	}
	return;
  }


  //
  // Calculates the dipole moment vector of a given water molecule according to the method of Morita/Hynes 2002
  //
  void MoritaH2O::SetDipoleMoment () {

	this->SetAtoms();	// first get the bonds set in the water

	// determine the bondlength displacements
	dR1 = this->OH1()->Magnitude() - BONDLENGTH_EQ;
	dR2 = this->OH2()->Magnitude() - BONDLENGTH_EQ;

	// angle displacement from the equilibrium
	dA = acos(this->Angle()) - ANGLE_EQ;

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
	  double distance = r.Magnitude();
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

  mpi::environment env(argc, argv);


  libconfig::Config cfg;
  cfg.readFile("system.cfg");

  std::string filename = cfg.lookup("analysis.morita2002.filename");
  libconfig::Setting &analysis = cfg.lookup("analysis");
  analysis.add("filename", libconfig::Setting::TypeString) = filename;

  WaterSystemParams wsp (&cfg, true);

  morita::SFGAnalyzer sfg (wsp);
  sfg.SystemAnalysis ();


  return 0;
}
