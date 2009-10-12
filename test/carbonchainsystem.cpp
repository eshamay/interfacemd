#include "../carbonchainsystem.h"

CarbonChainSystem::CarbonChainSystem
(
 const void * analysis_params,
 const WaterSystemParams& params
 )
: Analyzer(analysis_params, params)
{

	printf ("Loading up an Amber system with carbon chains to be analyzed:\n");
	// some info before starting
	printf ("The AMBER system loaded contains %d atoms\n", sys->size());

	return;
}

void CarbonChainSystem::Setup () {
  return;
}

void CarbonChainSystem::Analysis () {

  // find all the carbon-chain molecules
  FindMols (_ap->mol_name);

  // Here we need to find the histogram of the angles formed between the molecule's axis, and the reference axis for each timestep. 
  vector<int> timestep_histo = Molecular_Axis_Orientation_Histogram (
	  _ap->mol_name,
	  _ap->molecular_axis_function
	  )

	// update the running/total histogram
	RUN (timestep_histo) {
	  _ap->histogram[i] += timestep_histo[i];
	}

  return;
}

void CarbonChainSystem::DataOutput () {

  rewind (output);

  if (!(timestep % (output_freq * 10))) {
	// angle value is the row
	double angle;
	for (int ang = 0; ang < angbins; ang++) {
	  angle = ang * angres + angmin;
	  // print out the position for each position column
	  fprintf (output, "% 8.3f% 12d\n", angle, _ap->histogram[ang]);
	}
  }

  fflush (output);

  return;
}

void CarbonChainSystem::PostAnalysis () {
  // normalization of the histogram - account for the number of timesteps
  RUN (_ap->histogram) {
	_ap->histogram[i] /= timesteps;
  }

  return;
}

int main (int argc, const char **argv) {

  WaterSystemParams params;
  params.output = "C-chain.com.angle.dat";
  params.pbcflip = 30.0;
  params.output_freq = 500;

  struct AnalysisParams {
	string mol_name;
	AxisFunc molecular_axis_function;
	vector<int> histogram; 					// total running histogram
  }

  AnalysisParams ap;
  ap.mol_name = "dec";
  ap.molecular_axis_function = &CarbonChain::Vector_CoM_To_End;
  ap.histogram.resize(angbins, 0);

  CarbonChainSystem ccs ((void *)&ap, params);
  ccs.SystemAnalysis ();

  return 0;
}
