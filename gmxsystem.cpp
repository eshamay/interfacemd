#include "gmxsystem.h"

	namespace gromacs {
		GMXSystem::GMXSystem (const std::string gro_filepath, const std::string trr_filepath)
			:
				_gro(gro_filepath), _trr(trr_filepath)
		{ 
			this->_ParseMolecules();
		}

		void GMXSystem::LoadFirst() {
			_trr.LoadFirst();
			this->_ParseMolecules();
		}

		void GMXSystem::LoadNext() {
			_trr.LoadNext();
			this->_ParseMolecules();
		}

		void GMXSystem::_ParseMolecules () {

			VecR_it coord_it = _trr.begin_coords();
			//VecR_it vel_it = trr.begin_vels();
			VecR_it force_it = _trr.begin_forces();

			for (Atom_it it = _gro.begin(); it != _gro.end(); it++) {
				(*it)->Position(*coord_it);
				(*it)->Force(*force_it);

				++coord_it;
				//++vel_it;
				++force_it;
			}
		}

	} // namespace gromacs

/*
	 int main () {

	 GMXSystem sys ("sw_md.trr", "sw_md.gro");
	 Molecule * mol = sys.Molecules(0);
	 for (int i = 0; i < 10; i++) {
	 mol->Print();
	 sys.LoadNext();
	 }

	 return 0;
	 }
	 */
