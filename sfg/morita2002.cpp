#include "morita2002.h"

namespace morita {

  using namespace boost::numeric::ublas;

  SFGAnalyzer::SFGAnalyzer (WaterSystemParams& wsp)
	: 
	  Analyzer<AmberSystem> (wsp),
	  _p(0,0), _alpha(0,0)
  
  { return; }

  SFGAnalyzer::~SFGAnalyzer () {
	for (Morita_it it = _wats.begin(); it != _wats.end(); it++) { delete (*it); }
  }

  void SFGAnalyzer::Setup () {

	LoadWaters();
	int N = int_wats.size();

	_T.resize(N);
	_p.resize(N, 3);
	_alpha.resize(N, N);

  } // Setup

  void SFGAnalyzer::Analysis () {

	// Calculate the tensor 'T' which is formed of 3x3 matrix elements
	_T.clear();
	for (unsigned int i = 0; i < int_wats.size()/3; i++) {
	  for (unsigned int j = i+1; j < int_wats.size()/3; j++) {
		project(_T, slice(3*i,1,3), slice(3*j,1,3)) = math::DipoleFieldTensor(int_wats[i], int_wats[j]);
	  }
	}

	std::for_each (_wats.begin(), _wats.end(), SetDipoleMoment());
	_p.clear();
	for (unsigned i = 0; i < int_wats.size(); i++) {


	std::for_each (_wats.begin(), _wats.end(), SetPolarizability());


  } // Analysis

  void SFGAnalyzer::DataOutput (const unsigned int timestep) {

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
	_dipole.Zero();
	_dipole += _o->Position() * qO;
	_dipole += _h1->Position() * qH1;
	_dipole += _h2->Position() * qH2;
  }

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
	_DCM = this->DCMToLabMorita(z,1);
	_alpha = _alpha + (_DCM.Transpose() * _alpha1 * _DCM);

	_DCM = this->DCMToLabMorita(z,2);
	_alpha = _alpha + (_DCM.Transpose() * _alpha2 * _DCM);
  }


  


  namespace math {

	using namespace boost::numeric::ublas;

	template <typename T>
	void Tensor<T>::Print() const
	{ 
	  for (unsigned i = 0; i < this->size1(); i++) {
		for (unsigned j = 0; j < this->size2(); j++) {
		  printf ("% 8.3f", (*this)(i,j));
		}
		printf ("\n");
	  }
	  printf ("\n");
	}	// Print



	DipoleFieldTensor::DipoleFieldTensor (const Molecule* wat1, const Molecule* wat2)
	  :
		tensor_t(3,3)
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
