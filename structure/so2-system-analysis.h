#ifndef SO2_ANALYSIS_H_
#define SO2_ANALYSIS_H_

namespace so2_analysis {

	// a convenience class for working with systems comprised of at least 1 SO2 molecule and a whole bunch of waters
	template <typename T>
	class SO2SystemAnalyzer : public AnalysisSet< Analyzer<T> > {
		public:
			typedef Analyzer<T> system_t;

			SO2SystemAnalyzer (std::string desc, std::string filename) :
				AnalysisSet<system_t> (desc, filename) { }

			virtual ~SO2SystemAnalyzer () { 
				for (Wat_it it = this->all_wats.begin(); it != this->all_wats.end(); it++)
					delete *it;
				delete this->so2;
			}

			void Setup (system_t& t) {
				this->PreSetup(t);

				t.LoadAll();
				this->FindSO2 (t);
				this->so2->SetAtoms();
				this->s = so2->S();
				this->o1 = so2->O1();
				this->o2 = so2->O2();

				// grab the first location of the so2 as a reference for later analyses
				this->starting_point = system_t::Position(this->s);

				t.LoadWaters();
				// gather all the system waters into an analysis container
				for (Mol_it it = t.int_wats.begin(); it != t.int_wats.end(); it++) {
					WaterPtr wat (new Water(*(*it)));
					wat->SetAtoms();
					all_wats.push_back(wat);
				}

				all_wat_atoms.clear();
				std::copy(t.int_atoms.begin(), t.int_atoms.end(), std::back_inserter(all_wat_atoms));

				this->PostSetup(t);
			}

			// The method by which the SO2 of interest is found in the system
			virtual void FindSO2 (system_t& t) {
				// find the so2 in the system and set some pointers up
				MolPtr mol = Molecule::FindByType(t.sys_mols, Molecule::SO2);
				this->so2 = new SulfurDioxide(mol);
			}

			// things to do before the standard setup
			virtual void PreSetup (system_t& t) { }
			// things to do after the standard setup
			virtual void PostSetup (system_t& t) { }

			virtual void ReloadAnalysisWaters () {
				analysis_wats.clear();
				std::copy (all_wats.begin(), all_wats.end(), std::back_inserter(analysis_wats));
				// now all_wats has... all the waters, and analysis wats is used to perform some analysis
				analysis_atoms.clear();
				std::copy (all_wat_atoms.begin(), all_wat_atoms.end(), std::back_inserter(analysis_atoms));
			}


			// this assumes that the so2 is 'above' the surface. All waters in analysis wats above the so2 will be removed.
			virtual void FindWaterSurfaceLocation () {
				// get rid of everything above the so2
				this->analysis_wats.erase(
						remove_if(this->analysis_wats.begin(), this->analysis_wats.end(), system_t::MoleculeAbovePosition(this->starting_point, WaterSystem<T>::axis)), this->analysis_wats.end());

				// sort the waters by position along the reference axis - first waters are lowest, last are highest
				std::sort (this->analysis_wats.begin(), this->analysis_wats.end(), system_t::molecule_position_pred(Atom::O));
				int numWats = 20;				// number of waters to use for calculating the location of the "top" of the water surface
				surfaceLocation = 0.0;
				for (Wat_it it = this->analysis_wats.end() - 1; it != this->analysis_wats.end() - numWats; it--) {
					surfaceLocation += system_t::Position((*it)->ReferencePoint());
					//printf ("% .3f\n", system_t::Position((*it)->ReferencePoint()));
				}
				surfaceLocation /= numWats;
			}	// find surface water location


		protected:
			SulfurDioxide * so2;	// the sulfur dioxide of interest
			AtomPtr s, o1, o2;	// sulfur dioxide's atoms
			Water_ptr_vec all_wats, analysis_wats;
			Atom_ptr_vec all_wat_atoms, analysis_atoms;
			
			double starting_point;	// the original location of the so2 along the reference axis
			double surfaceLocation;	// location of the water surface along the reference axis
	};

}	// namespace so2_analysis
#endif
