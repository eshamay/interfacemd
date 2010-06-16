#include "test.h"

Test::Test (WaterSystemParams& wsp)
  :	Analyzer<AmberSystem> (wsp)
{

  this->sys = new AmberSystem("prmtop", "mdcrd", "mdvel");

  return;
}

void Test::Setup () {

  //this->LoadAll();
  //this->LoadWaters();


  //std::for_each(int_atoms.begin(), int_atoms.end(), tform);

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

  std::vector<int> v;
  for (int i = 0; i < 10; i++) {
	v.push_back(i*10);
  }

  for (std::vector<int>::const_iterator it = v.begin(); it != v.end(); it++)
	cout << *it << " ";
  cout << endl;

  v.erase(std::remove_if(v.begin(), v.end(), std::not1(std::bind1st(std::less<int>(), 30))), v.end());
  
  for (std::vector<int>::const_iterator it = v.begin(); it != v.end(); it++)
	cout << *it << " ";
  cout << endl;

  //t.SystemAnalysis ();

  return 0;
}

