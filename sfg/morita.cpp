#include "morita.h"

namespace morita {

	SFGCalculator MoritaAnalysis::sfg;
	Complex_vec MoritaAnalysis::TimestepChi;
	Complex_vec MoritaAnalysis::TotalChi;

	unsigned long MoritaAnalysis::numMolsProcessed = 0;
	bool MoritaAnalysis::firstMol = true;
	bool MoritaAnalysis::firstTimeStep = true;

	MoritaAnalysis::~MoritaAnalysis () {
		for (Morita_it it = analysis_wats.begin(); it != analysis_wats.end(); it++) {
			delete *it;
		}
	}

	void MoritaAnalysis::SetupSystemWaters (system_t& t) {

		// load all the waters into the int_wats container
		t.LoadWaters();
		//std::for_each(t.int_wats.begin(), t.int_wats.end(), std::mem_fun(&Molecule::Print));

		for (Morita_it it = all_wats.begin(); it != all_wats.end(); it++) {
			delete *it;
		}
		all_wats.clear();

		// load up the all_wats with new derived Morita waters that have some extra functionality
		for (Mol_it it = t.int_wats.begin(); it != t.int_wats.end(); it++) {
			MoritaH2O_ptr ptr (new MoritaH2O (*it));
			all_wats.push_back(ptr);
		}

		// analysis_wats will hold only those waters that will by analyzed
		analysis_wats.clear();
		std::copy(all_wats.begin(), all_wats.end(), back_inserter(analysis_wats));
		// here filter out the waters to use for analysis
		this->SelectAnalysisWaters ();

		return;
	}

	void MoritaAnalysis::SelectAnalysisWaters () {
		// sort by position in the slab
		std::sort(analysis_wats.begin(), analysis_wats.end(), Analyzer<AmberSystem>::molecule_position_pred(Atom::O));

		// then slice the waters to keep only a certain number of them
		int numAnalysisWaters = WaterSystem<AmberSystem>::SystemParameterLookup ("analysis.morita2002.number-of-analysis-waters");
		analysis_wats.erase (analysis_wats.begin(), analysis_wats.end() - numAnalysisWaters);
	}

	void MoritaAnalysis::Analysis (system_t& t) {

		printf ("test\n");
		TimestepChi.clear();	// it's a new timestep
		firstMol = true;		// every timestep we will have to go through all the molecules again

		this->SetupSystemWaters (t);

		// and then update our bond data to reflect the interfacial region and find all the hydrogen bonds
		t.UpdateGraph ();

		std::for_each(analysis_wats.begin(), analysis_wats.end(), this->SFGProcess);

		// we collect the data for each timestep into the running total for the end spectrum
		CollectChi (TotalChi.begin(), TotalChi.end(), TimestepChi.begin());

		return;
	}

	// output data to the file
	void MoritaAnalysis::DataOutput (system_t& t) {

		rewind (t.Output());

		double freq = START_FREQ * sfg_units::AU2WAVENUMBER;
		double scale = (double)numMolsProcessed;
		double r, im;

		for (Complex_vec::const_iterator it = TotalChi.begin(); it != TotalChi.end(); it++) {
			r = real(*it) / scale;
			im = imag(*it) / scale;

			fprintf (t.Output(), "% 20.4e% 20.4e% 20.4e\n", freq, r, im);

			freq += FREQ_STEP * sfg_units::AU2WAVENUMBER;
		}

		fflush (t.Output());

		return;
	}

} // namespace morita


int main () {

	Analyzer<AmberSystem> analyzer;
	morita::MoritaAnalysis analysis;
	analyzer.SystemAnalysis(analysis);

	return 0;
}

