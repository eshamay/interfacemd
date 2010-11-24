#ifndef __RDF_ANALYSIS_H
#define __RDF_ANALYSIS_H

#include "../analysis.h"

namespace md_analysis {

	template <class T>
		class rdf_analysis : public AnalysisSet< Analyzer<T> > {
			public:
			typedef Analyzer<T> system_t;
				rdf_analysis () :
					AnalysisSet<system_t> (
						std::string("RDF Analysis"),
						std::string("rdf.dat")) { }

				virtual ~rdf_analysis () {
					//delete so2;
				}
				void Setup (system_t& t);
				void Analysis (system_t& t);
				void DataOutput (system_t& t);

			protected:
				typedef std::pair<Atom::Element_t, Atom::Element_t>	element_pair_t;
				// a map of atom types/elements pairs to histograms - 1 histogram for each atom type pair to be analyzed
				typedef std::map<element_pair_t, histogram_utilities::Histogram1D<double> >	element_histogram_map_t;
				element_histogram_map_t histograms;
				// holds the list of element pairs to be analyzed
				typedef std::list<element_pair_t> element_pair_list;
				element_pair_list	element_pairs;
				//SulfurDioxide * so2;
		};

	template <typename T>
	void rdf_analysis<T>::Setup(system_t& t) {
		// get the rdf parameters from the configuration file
		double minimum = t.SystemParameterLookup("analysis.rdf.minimum");
		double maximum = t.SystemParameterLookup("analysis.rdf.maximum");
		double resolution = t.SystemParameterLookup("analysis.rdf.resolution");

		// form all the histograms and the mapping between element pairs with the histograms to be used for rdf binning
		libconfig::Setting &atompairs = t.SystemParameterLookup("analysis.rdf.atom-pairs");
		for (int i = 0; i < atompairs.getLength(); i++) {
			// grab all the atom name-pairs from the configuration file and form the list of element pairs to be RDF'ed
			element_pair_t ep = std::make_pair(Atom::String2Element(atompairs[i][0]), Atom::String2Element(atompairs[i][1]));
			element_pairs.push_back (ep);
			// secondly insert the new mapping into our collection
			histograms.insert (std::make_pair(ep, histogram_utilities::Histogram1D<double>(minimum, maximum, resolution)));
		}

		// load all the atoms to work with for the duration of the analysis
		t.LoadAll();

		/*
		// grab the so2 of interest - for now!
		int id = t.SystemParameterLookup ("analysis.reference-molecule-id");
		MolPtr mol = Molecule::FindByID (t.sys_mols, id);
		so2 = new SulfurDioxide(mol);
		so2->SetAtoms();
		*/
	} // rdf analysis setup

	template <typename T>
	void rdf_analysis<T>::Analysis(system_t& t) {
		double atomic_distance;
		// go through each atom pair and 
		for (Atom_it it = t.sys_atoms.begin(); it != t.sys_atoms.end() - 1; it++) {
			for (Atom_it jt = it+1; jt != t.sys_atoms.end(); jt++) {
				// first check if the pair of atoms is one of the pairs to analyze - is the pair in the list of pairs?
				element_pair_t ep = std::make_pair((*it)->Element(), (*jt)->Element());
				element_pair_list::iterator epl_it = pair_utility::PairListMember (ep, element_pairs.begin(), element_pairs.end());

				if (epl_it == element_pairs.end()) continue;	// if the pair isn't in the list, then don't bin it!

				// for each pair in the list, bin the distance between the atoms into the proper histogram
				atomic_distance = MDSystem::Distance(*it, *jt).norm();
				histograms.find(ep)->second(atomic_distance);
			}
		}

		/*
		// specific analysis to bag only the O1 and O2 rdfs for so2
		t.LoadWaters();
		Atom::KeepByElement (t.int_atoms, Atom::H);
		AtomPtr o1 = so2->O1();
		AtomPtr o2 = so2->O2();
		for (Atom_it it = t.int_atoms.begin(); it != t.int_atoms.end(); it++) {
				atomic_distance = MDSystem::Distance(*it, o1).norm();
				histograms.begin()->second(atomic_distance);
				atomic_distance = MDSystem::Distance(*it, o2).norm();
				histograms.begin()->second(atomic_distance);
		}
		*/

	}	// rdf analysis analysis



	template <typename T>
		void rdf_analysis<T>::DataOutput(system_t& t) {
			rewind (t.Output());

			double minimum = t.SystemParameterLookup("analysis.rdf.minimum");
			double maximum = t.SystemParameterLookup("analysis.rdf.maximum");
			double resolution = t.SystemParameterLookup("analysis.rdf.resolution");

			double dV, n, N;
			double total_volume = 4.0/3.0 * M_PI * pow(maximum, 3);

			// go through each value in the RDF position range
			for (double r = minimum; r < maximum; r += resolution) {
				// print the position to the first column
				fprintf (t.Output(), "%12.3f ", r);
				dV = 4.0 * M_PI * pow(r, 2) * resolution;	// volume of the shell being considered

				// then for each element pair analyzed, print a column with the value of g(r) for the given position
				for (element_histogram_map_t::iterator it = histograms.begin(); it != histograms.end(); it++) {
					n = it->second.Population(r);	// the population of the shell
					N = it->second.Count();				// total population of the histogram

					fprintf (t.Output(), "%12.3f ", n * total_volume / dV / N);
				}
				fprintf(t.Output(), "\n");
			}

			fflush(t.Output());
			return;
		}	// rdf analysis data output

} // namespace md_analysis
#endif
