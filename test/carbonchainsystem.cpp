#include "carbonchainsystem.h"

	CarbonChainSystem::CarbonChainSystem (int argc, const char **argv, const WaterSystemParams& params)
:	WaterSystem<AmberSystem>(params)
	
{

	printf ("Loading up an Amber system with carbon chains to be analyzed:\n");

	this->sys = new AmberSystem("prmtop", "mdcrd", "mdvel");

	// some info before starting
	printf ("The AMBER system loaded contains %d atoms\n", sys->size());

	return;
}

// Create a histogram of the angles formed by an axis of a carbon-chain and a given reference axis. The supplied axis function (pointer) should return the molecular axis vector of the carbon-chain.
vector<int> CarbonChainSystem::OrientationHistogram (
		const WaterSystemParams& params,
		const vector<Molecule *> mols,
		const string name,
		VecR (*axisFunc)(const Atom * atom))
{

	// set up the histogram for output
	vector<int> histo (anglebins, 0);

	// Run an analysis on all the carbon-chain molecules in the system to find their orientations over the course of a simulation with respect to a given axis.
	RUN (mols) {
		// Search for the named molecules in the system
		if (mols[i]->Name().find(name) == string::npos) continue;

		// find the particular molecular-axis vector
		Atom * c10 = mols[i]->Carbon(9);
		VecR molAxis = mols[i]->axisFunc(c10);
		double angle = (molAxis < params.ref_axis);

		int anglebin = int ((angle - params.anglemin)/params.angleres);
		histo[anglebin]++;
	}

	return (histo);
}


void CarbonChainSystem::Orientation-Analysis () {

	string mol_name = "dec";
	typedef mol_t Decane;
	typedef vector<int>(*axisFunc)(const Atom *);
	axisFunc = &CarbonChain::Vector_CoM_To_Atom;

	// total running histogram
	vector<int> histogram (angbins, 0);

	// start the analysis - run through each timestep
	for (timestep = 0; timestep < timesteps; timestep++) {

		// find all the carbon-chain molecules
		this->FindMols (mol_name);

		// Here we need to find the histogram for each timestep of the angles formed between the molecule's axis, and the reference axis
		vector<int> timestep_histo = 
			dec->OrientationHistogram(
					_params,
					int_mols,
					mol_name,
					axisFunc)

			// update the running/total histogram
			RUN (timestep_histo) {
				histogram[i] += timestep_histo[i];
			}

		this->sys->LoadNext();

		this->OutputStatus ();
		this->Output_Orientation_Histogram_Data (histogram);
	}

	RUN (histogram) {
		histogram[i] /= timesteps;
	}

	this->Output_Orientation_Histogram_Data (histogram);

	return;
}




int main (int argc, const char **argv) {

	WaterSystemParams params;
	params.output = "decane.C-chain.angle.dat";
	params.pbcflip = 30.0;
	params.output_freq = 500;

	CarbonChainSystem ccs (argc, argv, params);

	ccs.Orientation-Analysis();

	return 0;
}
