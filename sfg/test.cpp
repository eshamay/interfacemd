#include "test.h"

Test::Test (WaterSystemParams& wsp)
  :	Analyzer<AmberSystem> (wsp)
{

  this->sys = new AmberSystem("prmtop", "mdcrd", "mdvel");

  return;
}

void Test::Setup () {

  this->LoadAll();
  this->LoadWaters();

  class tform : public std::unary_function<Atom *,Atom *> {
	public:
	  Atom * operator() (Atom * atom) {
		atom->Print();
		return atom;
	  }
  };

  std::vector<Atom *> bools;
  std::transform(
	  int_atoms.begin(), 
	  int_atoms.end(), 
	  bools.begin(), 
	  tform()
	  );

  //Water * wat = static_cast<Water *>(int_wats[0]);
  //wat->Print();
  return;
}

int main () {

  libconfig::Config cfg;
  cfg.readFile("system.cfg");

  std::string filename = cfg.lookup("analysis.test.filename");
  libconfig::Setting &analysis = cfg.lookup("analysis");
  analysis.add("filename", libconfig::Setting::TypeString) = filename;

  WaterSystemParams wsp (cfg);

  Test t (wsp);

  //t.SystemAnalysis ();

  return 0;
}

