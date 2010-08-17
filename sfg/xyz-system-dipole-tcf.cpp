#include "xyz-analyzer.h"

void XYZAnalyzer::Analysis () {

  LoadAll();

  // calculate the dipole for each molecule
  std::for_each (sys_mols.begin(), sys_mols.end(), MDSystem::CalcDipole);

  // grab all the dipoles into a single container
  VecR_vec dipoles;
  std::transform (sys_mols.begin(), sys_mols.end(), std::back_inserter(dipoles), std::mem_fun<VecR,Molecule>(&Molecule::Dipole));

  // Sum all the dipoles for the total system dipole
  VecR M = std::accumulate (dipoles.begin(), dipoles.end(), VecR());

  if (!timezero) {
	M_0 = M;
	timezero = true;
  }

  fprintf (this->output, "%f %f %f\n", M(0), M(1), M(2));
  //fprintf (this->output, "% 20.6f\n", M_0.dot(M));
  //for (Mol_it mol = sys_mols.begin(); mol != sys_mols.end(); mol++) {
	//M += this->sys->CalcDipole(*mol);
  //}

  return;
}	// Analysis

  /*
void XYZAnalyzer::DataOutput () {

  // process time-domain data
  std::vector<double> corr;
  VecR M_0 = *_M.begin();
  for (VecR_it it = _M.begin(); it != _M.end(); it++)
	corr.push_back (M_0.dot(*it));

  rewind(output);
  // print the time-domain data and also the fft version
  for (std::vector<double>::const_iterator it = corr.begin(); it != corr.end(); it++) {
	fprintf (output, "% 20.6f\n", *it);
  }
  */

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

  return; 
}	// data output
   */



int main () {

  libconfig::Config cfg;
  cfg.readFile("system.cfg");

  std::string filename ("xyz-system-dipole-vector.dat");
  libconfig::Setting &analysis = cfg.lookup("analysis");
  analysis.add("filename", libconfig::Setting::TypeString) = filename;

  WaterSystemParams wsp (cfg);

  XYZAnalyzer test (wsp);

  test.SystemAnalysis ();

  return 0;
}

