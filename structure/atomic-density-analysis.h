#ifndef __DENSITY_HISTOGRAM_H
#define __DENSITY_HISTOGRAM_H

#include "../utility.h"
#include "../analysis.h"

/*
 * This analysis package calculates atomic densities throughout a system along the 3 different axes of the system (X,Y,Z). A system.cfg file is required in the working directory, and the input files from the system must have specific names (or links to the files with specific names). For Amber Systems there must be a files.prmtop file (or link) and a files.mdcrd file, using those names exactly as lsited in system.cfg.
 *
 * In the system.cfg file, aside from all the other usual (mandatory) sections, another section is required for this analysis:
 *
 * analysis.density.atom-names
 *		This section is a list of the atom names (as defined in the xyz or prmtop files) that will be analyzed and output to the data file.
 *		example: 
 *				density:
					{
						atom-names = ( "O", "H1", "H2" );
					};

 * Output will be written to the file specific in system.cfg's analysis.density.filename and will overwrite previously written data if not backed up.
 *
 * Data in the file will be a series of columns, 3 column-pairs for x,y,z directions for each atom listed. Each column pair contains the position along the direction, and the number density of the atoms. Thus, for the example listed above, the columns would be:
 *		X-pos, O_x_population, Y-pos, O_y_population, Z-pos, O_z_population, X-pos, H1_x_population, etc.
 *
 *		Analysis will use the position extents listed in system.cfg's analysis.position-range +/- 10-angstroms, with a slice resolution of analysis.resolution.position.
 */

namespace md_analysis {

	template <typename T>
		class atomic_density_analysis : public AnalysisSet< Analyzer<T> > {
			public:
				typedef Analyzer<T> system_t;

				atomic_density_analysis () : 
					AnalysisSet<system_t> (
							std::string("An analysis of the density of atoms in a system based on atomic position"),
							//std::string("3d-atomic-density.dat")) { }
							WaterSystem<T>::SystemParameterLookup("analysis.density.filename")) { }

				virtual ~atomic_density_analysis () { }

				void Setup (system_t& t);
				void Analysis (system_t& t);
				// For each atom type (name) in the system, the histograms in each direction will be output
				void DataOutput (system_t& t);

			protected:
				std::vector<std::string> atom_name_list;
				// Every atom-name will have its own histogram of positions in the system. Each position is held as a vector to the atom site.
				typedef histogram_utilities::Histogram1D<double>	histogram_t;
				typedef std::vector<histogram_t>									histogram_set;
				typedef std::pair<std::string, histogram_set>			histogram_map_elmt;
				typedef std::map<std::string, histogram_set>			histogram_map;
				histogram_map histograms;

				class atomic_position_binner : public std::binary_function<AtomPtr,histogram_map *,void> {
					public:
						void operator() (AtomPtr atom, histogram_map * histos) const {
							histogram_set& hs = histos->operator[] (atom->Name());
							hs[0].operator() (atom->X());
							hs[1].operator() (atom->Y());
							hs[2].operator() (atom->Z());
						}
				} binner;

		};


	template <class T>
		void atomic_density_analysis<T>::Setup (system_t& t) {

			AnalysisSet<system_t>::Setup(t);

			// grab the list of atomic names/types that will be used for the analysis and create the vector-histograms
			libconfig::Setting &atom_names = WaterSystem<T>::SystemParameterLookup("analysis.density.atom-names");
			for (int i = 0; i < atom_names.getLength(); i++)
			{
				std::string atom_name = atom_names[i];
				atom_name_list.push_back(atom_name);

				histogram_set hs (3, histogram_t(-10.0, WaterSystem<T>::posmax+10.0, WaterSystem<T>::SystemParameterLookup("analysis.resolution.position")));
				histograms.insert(histogram_map_elmt(atom_name, hs));
			}

			// narrow down the system atoms to just those with names we're looking for
			md_name_utilities::KeepByNames (t.int_atoms, atom_name_list);

		}	// Setup


	template <class T>
		void atomic_density_analysis<T>::Analysis (system_t& t) { 

			//for (Atom_it it = t.int_atoms.begin(); it != t.int_atoms.end(); it++) {
			//VecR pos = (*it)->Position();
			//histograms[(*it)->Name()].push_back (pos);
			//}

			std::for_each (t.int_atoms.begin(), t.int_atoms.end(), std::bind2nd(binner, &histograms));

		}

	template <class T>
		void atomic_density_analysis<T>::DataOutput (system_t& t) {

			rewind(t.Output());

			// first output the header of all the atom-names
			for (std::vector<std::string>::const_iterator it = atom_name_list.begin(); it != atom_name_list.end(); it++) {
				fprintf (t.Output(), "%s_x %s_y %s_z ", it->c_str(), it->c_str(), it->c_str());
			}
			fprintf (t.Output(), "\n");
			fflush(t.Output());

			// output the data from the histograms
			double dr = histograms.begin()->second[0].Resolution();
			double min = histograms.begin()->second[0].Min();
			double max = histograms.begin()->second[0].Max();

			for (double r = min; r < max; r+=dr) {
				for (std::vector<std::string>::const_iterator name = atom_name_list.begin(); name != atom_name_list.end(); name++) {

					histogram_set& hs = histograms[*name];
					for (int ax = 0; ax < 3; ax++) {
						fprintf (t.Output(), "% 8.3f % 8.3f ", r, hs[ax].Population(r)/t.timestep);
					}
				}
				fprintf (t.Output(), "\n");
			}

			fflush(t.Output());

		}	// Data Output


	/*
	 * This analysis is for an SO2/H2O system where the SO2 starts outside of the water and then adsorbs to the surface and is absorbed into the water phase. To watch the adsorption rate, the water phase is analyzed and the topmost and lowest water (along the main axis) is determined and then all the SO2 molecules located within the water phase (and a fudge-factor if desired) are counted up and printed out.
	 * The output will be 1 column in a datafile:
	 *		# of SO2s in the water phase
	 */
	template <typename T>
		class so2_uptake_analysis : public AnalysisSet< Analyzer<T> > {
			public:
				typedef Analyzer<T> system_t;

				so2_uptake_analysis () : 
					AnalysisSet<system_t> (
							std::string("An analysis of the adsorption of SO2 into a water phase"),
							//std::string("3d-atomic-density.dat")) { }
							std::string("so2-adsorption.dat")) { }

				virtual ~so2_uptake_analysis () { }

				void Analysis (system_t& t);
			protected:
				MolPtr high_water, low_water;
				double high_position, low_position;
				int numAdsorbed;

				void FindHighAndLowWaters (system_t& t);
				void CountSO2InWater (system_t& t);
		};

	template <typename T>
	void so2_uptake_analysis<T>::FindHighAndLowWaters (system_t& t) {

		t.LoadWaters();

		// sort the waters according to position in the slab
		std::sort(t.int_wats.begin(), t.int_wats.end(), Analyzer<AmberSystem>::molecule_position_pred(Atom::O));
		// find the highest and the lowest at the extents of the water slab
		high_water = t.int_wats.back();
		high_position = Analyzer<T>::Position(high_water->GetAtom(Atom::O));
		low_water = t.int_wats.front();
		low_position = Analyzer<T>::Position(low_water->GetAtom(Atom::O));

	}

	template <typename T>
	void so2_uptake_analysis<T>::CountSO2InWater (system_t& t) {

		// get the listing of the SO2s in the system
		t.LoadAll();
		t.int_mols.clear();
		for (Mol_it so2 = t.sys_mols.begin(); so2 != t.sys_mols.end(); so2++) {
			if ((*so2)->MolType() == Molecule::SO2)
				t.int_mols.push_back(*so2);
		}

		Atom * s;
		numAdsorbed = 0;
		double pos;
		for (Mol_it so2 = t.int_mols.begin(); so2 != t.int_mols.end(); so2++) {
			pos = (*so2)->GetAtom(Atom::O)->Position()[WaterSystem<T>::axis];
			if (pos < high_position && pos > low_position) {
				++numAdsorbed;
			}
		}
	}

	template <typename T>
	void so2_uptake_analysis<T>::Analysis (system_t& t) {
		this->FindHighAndLowWaters(t);
		this->CountSO2InWater(t);

		fprintf (t.Output(), "%12d %8d\n", t.Timestep(), this->numAdsorbed);
	}

}	// namespace md_analysis



#endif
