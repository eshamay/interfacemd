#include "test.h"

Tester::Tester (WaterSystemParams& wsp)
  :	
	Analyzer<AmberSystem> (wsp)
{

  return;
}

void Tester::Setup () {
  LoadWaters ();

  WaterPtr wat (new Water (int_wats[0]));

  wat->Print();
  wat->SetAtoms();
  VecR oh1 = *wat->OH1();
  VecR oh2 = *wat->OH2();

  printf ("oh bond 1\n");
  oh1.Print();
  MatR dcm = wat->DCMToLabMorita(z,1);

  MatR Alab;
  Alab.Set(0,0,1.0);
  Alab.Set(1,2,4.0);
  Alab.Set(2,1,0.3);

  VecR Blab = Alab * oh1;
  printf ("B lab\n");
  Blab.Print();

  //MatR Aloc = dcm.Transpose() * Alab * dcm;
  MatR Aloc = dcm * Alab * dcm.Transpose();
  VecR oh1loc = dcm.Transpose() * oh1;
  printf ("oh1 loc\n");
  oh1loc.Print();

  VecR Bloc = Aloc * oh1loc;
  printf ("B loc\n");
  Bloc.Print();

  printf ("A*b in the loc then rotated back\n");
  (dcm * Bloc).Print();


  return;
}

void Tester::Analysis () {

  return;
}

void Tester::DataOutput (const unsigned int timestep) { 

  return; 
}

int main () {

  libconfig::Config cfg;
  cfg.readFile("system.cfg");

  std::string filename = cfg.lookup("analysis.test.filename");
  libconfig::Setting &analysis = cfg.lookup("analysis");
  analysis.add("filename", libconfig::Setting::TypeString) = filename;

  WaterSystemParams wsp (cfg);

  Tester test (wsp);

  test.SystemAnalysis ();

  return 0;
}

