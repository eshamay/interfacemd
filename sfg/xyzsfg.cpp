#include "xyzsfg.h"


XYZSFGAnalyzer::XYZSFGAnalyzer (WaterSystemParams& wsp)
  :	
	Analyzer<XYZSystem> (wsp)
{
  return;
}

void XYZSFGAnalyzer::Setup () {

  return;
}

void XYZSFGAnalyzer::Analysis () {

  LoadAll();

  VecR M (0.0,0.0,0.0);

  for (Mol_it mol = sys_mols.begin(); mol != sys_mols.end(); mol++) {
	M += this->sys->CalcDipole(*mol);
  }
  _M.push_back(M);

  return;
}

void XYZSFGAnalyzer::DataOutput (const unsigned int timestep) { 

  // process time-domain data
  std::vector<double> corr;
  for (VecR_it it = _M.begin(); it != _M.end(); it++)
	corr.push_back (_M[0] * *it);

  // print the time-domain data and also the fft version
  for (std::vector<double>::const_iterator it = corr.begin(); it != corr.end(); it++) {
	fprintf (output, "% 20.6f\n", *it);
  }

  // ********  fftw work ******** //
  /*
	 fftw_complex *in;
	 fftw_complex *corr_fft;
	 fftw_plan p;

	 in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*corr.size());

  // sps transform
  for (int i = 0; i < corr.size(); i++) {
  in[i][0] = corr[i];
  in[i][1] = 0.0;
  }
  //std::copy(corr.begin(), corr.end(), in);
  corr_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*corr.size());
  p = fftw_plan_dft_1d(corr.size(),in,corr_fft,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_execute(p);

  // output one column each for
  // time-domain data, real, imaginary, magnitude
  rewind(output);
  for (int i = 0; i < corr.size(); i++) {
  fprintf (output, "% 20.6f % 20.6f % 20.6f % 20.6f\n", 
  corr[i], 
  corr_fft[i][0], 
  corr_fft[i][1], 
  sqrt(corr_fft[i][0]*corr_fft[i][0] + corr_fft[i][1]*corr_fft[i][1]));
  }

  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(corr_fft);
   */

  return; 
}



int main () {

  libconfig::Config cfg;
  cfg.readFile("system.cfg");

  std::string filename = cfg.lookup("analysis.xyzsfg.filename");
  libconfig::Setting &analysis = cfg.lookup("analysis");
  analysis.add("filename", libconfig::Setting::TypeString) = filename;

  WaterSystemParams wsp (cfg);

  XYZSFGAnalyzer test (wsp);

  test.SystemAnalysis ();

  return 0;
}

