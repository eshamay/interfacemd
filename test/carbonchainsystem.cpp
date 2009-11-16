#include "carbonchainsystem.h"

CarbonChainSystem::CarbonChainSystem
(
 AnalysisParams& analysis_params,
 WaterSystemParams& params
 )
: Analyzer<CarbonChainSystem, AnalysisParams> (analysis_params, params)
{

	printf ("Loading up an Amber system with carbon chains to be analyzed:\n");
	// some info before starting
	printf ("The AMBER system loaded contains %d atoms\n", sys->size());

	return;
}

void CarbonChainSystem::EmptyFn () {
	return;
}

void CarbonChainSystem::ClearHistogram_FindMols () {
	FindMols (_ap.mol_name);

	RUN (_ap.density_histo_2) 
		_ap.density_histo_2[i].clear();
	_ap.density_histo_2.clear();

	RUN (_ap.histogram_2) 
		_ap.histogram_2[i].clear();
	_ap.histogram_2.clear();

	return;
}

void CarbonChainSystem::MolecularAxisDirection () {

	// Here we need to find the histogram of the angles formed between the molecule's axis, and the reference axis for each timestep. 
	vector<int> timestep_histo = 
		Molecular_Axis_Orientation_Histogram<mol_t> (
				_ap.mol_name,
				_ap.axisFn
				);

	// update the running/total histogram
	RUN (timestep_histo) {
		_ap.histogram[i] += timestep_histo[i];
	}

	return;
}

void CarbonChainSystem::AngleHistogramOutput (const int timestep) {

	rewind (output);

	if (!(timestep % (output_freq * 10))) {
		// angle value is the row
		double angle;
		for (int ang = 0; ang < angbins; ang++) {
			angle = ang * angres + angmin;
			// print out the position for each position column
			fprintf (output, "% 8.3f% 12.3f\n", angle, _ap.histogram[ang]);
		}
	}

	fflush (output);

	return;
}

void CarbonChainSystem::Histogram_2_Output (const int timestep) {

	rewind (output);

	if (!(timestep % (output_freq * 10))) {
		double val = 0.0;
		RUN(_ap.histogram_2) {
			RUN2 (_ap.histogram_2[i]) {
				if (_ap.density_histo_2[i][j] == 0)
					val = 0.0;
				else {
					val = _ap.histogram_2[i][j]/(double)_ap.density_histo_2[i][j];
				}
				fprintf (output, "% 8.3f", val);
			}
			fprintf (output, "\n");
		}
	}

	fflush (output);

	return;
}

void CarbonChainSystem::Histogram_2_Output_2 (const int timestep) {

	rewind (output);

	if (!(timestep % (output_freq * 10))) {
		double val = 0.0;
		RUN(_ap.density_histo_2) {
			RUN2 (_ap.density_histo_2[i]) {
				val = (double)_ap.density_histo_2[i][j]/(double)timestep;
				fprintf (output, "% 9.4f", val);
			}
			fprintf (output, "\n");
		}
	}

	fflush (output);



	return;
}

void CarbonChainSystem::MolecularPlane () {

	vector< vector<double> > timestep_histo = 
		Interface_Location_Histogram<mol_t> (
				x, z, 
				_ap.mol_name, 
				_ap.atom_name,
				_ap.histogram_2,
				_ap.density_histo_2);

	return;
}

void CarbonChainSystem::HistogramAveraging () {
	// normalization of the histogram - account for the number of timesteps
	RUN (_ap.histogram) {
		_ap.histogram[i] /= (double)timesteps;
	}

	return;
}

void CarbonChainSystem::HistogramAveraging_2 () {
	// normalization of the histogram - account for the number of timesteps
	RUN (_ap.histogram_2) {
		RUN2 (_ap.histogram_2) {
			_ap.histogram_2[i][j] /= (double)timesteps;
		}
	}

	return;
}

int main (int argc, const char **argv) {

	WaterSystemParams params ("C10-locationdensity.dat");
	params.pbcflip = 30.0;
	params.output_freq = 500;

	AnalysisParams ap;
	ap.mol_name = "pds";
	ap.atom_name = "C10";
	ap.histogram.resize(params.angbins+1, 0.0);
	ap.axisFn = &CarbonChain::Vector_CoM_To_End;

	CarbonChainSystem ccs (ap, params);
	ccs.Set_Setup (&CarbonChainSystem::ClearHistogram_FindMols);
	ccs.Set_Analysis (&CarbonChainSystem::MolecularPlane);
	ccs.Set_PostAnalysis (&CarbonChainSystem::EmptyFn); 
	ccs.Set_DataOutput (&CarbonChainSystem::Histogram_2_Output_2);

	ccs.SystemAnalysis (&ccs);

	return 0;
}
