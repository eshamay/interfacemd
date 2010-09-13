#include "amber-morita2002.h"

namespace morita {

	void AmberMorita2008Analysis::SelectAnalysisWaters () {

		// slice the waters in the system according to the center-of-mass location
		/*
		// find the system center of mass
		VecR com = Analyzer<AmberSystem>::CenterOfMass (analysis_wats.begin(), analysis_wats.end());
		double cutoff = com[WaterSystem<AmberSystem>::axis] + 20.0;

		std::cout << "cutoff set to: " << cutoff << std::endl;
		*/

		// sort the waters in the system according to position along the main axis
		std::sort(analysis_wats.begin(), analysis_wats.end(), Analyzer<AmberSystem>::molecule_position_pred(Atom::O));

		/*
		for (Morita_it it = analysis_wats.begin(); it != analysis_wats.end(); it++) {
			(*it)->GetAtom(Atom::O)->Print();
		}
		*/


		// only use waters found above a particular location in the slab (i.e. center of mass)
		//std::pair<double,double> slice = std::make_pair(cutoff,WaterSystem<AmberSystem>::posmax);
		//WaterSystem<AmberSystem>::SliceWaters<MoritaH2O_ptr> (analysis_wats, slice);

		//std::cout << analysis_wats.size() << " waters to be analyzed" << std::endl;
		analysis_wats.erase (analysis_wats.begin(), analysis_wats.end() - 200);
		//std::cout << analysis_wats.size() << " waters to be analyzed" << std::endl;

		return;
	}

	void AmberMorita2008Analysis::SetAnalysisWaterDipoleMoments () {
		std::for_each (analysis_wats.begin(), analysis_wats.end(), std::mem_fun(&MoritaH2O::SetDipoleMoment));
		return;
	}

	void AmberMorita2008Analysis::SetAnalysisWaterPolarizability () {

		MatR alpha;
		MatR dcm;

		// for every water to be analyzed:
		//		set the molecule up
		//		calculate the direction cosine matrix from the local to lab frame rotation
		//		find the polarizability from the lookup-table
		//		rotate the polarizability to the lab frame
		//		set the molecule's polarizability to the value
		for (Morita_it it = analysis_wats.begin(); it != analysis_wats.end(); it++) {

			(*it)->SetAtoms();	// first get the atoms and bonds set in the water

			dcm = MatR::Zero();
			dcm = (*it)->DCMToLabMorita(1);	// set up the direction cosine matrix for rotating the polarizability to the lab-frame

			// lookup the polarizability from the data file
			alpha = MatR::Zero();
			alpha = pdf.Matrix(
						(*it)->OH1()->Magnitude(), 
						(*it)->OH2()->Magnitude(), 
						acos((*it)->Angle()) * 180.0/M_PI);
			
			// rotate the polarizability tensor into the lab-frame
			alpha = dcm * alpha;

			(*it)->SetPolarizability (alpha);
		}

	}	// set polarizability

} // namespace morita





// used to calculate the SFG spectrum based on the morita/hynes 2008 method
int main (int argc, char **argv) {

	Analyzer<AmberSystem> analyzer;
	morita::AmberMorita2008Analysis analysis;
	analyzer.SystemAnalysis(analysis);

	return 0;
}

