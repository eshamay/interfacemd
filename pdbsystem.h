#ifndef PDBFILE_H_
#define PDBFILE_H_

#include "mdfiles.h"
#include "moleculefactory.h"
#include <cstring>
#include <boost/algorithm/string.hpp>


namespace md_files {

	class PDBSystem : public MDSystem, public CoordinateFile {

		protected: 
			int _numAtoms;
			bool _initialized;

			void _ParseAtoms ();
			void _ParseMolecules ();

			void CreateAtoms();
			void CountAtoms ();

		public:

			PDBSystem (std::string path) : CoordinateFile(path) {
				LoadFirst();
			}

			void LoadFirst () {
				if (!_initialized) {
					CountAtoms();
					CreateAtoms();
					_initialized = true;
				}

				rewind (_file);
				LoadNext();
			}

			// Various control functions
			void LoadNext ();

			//void WritePDB ();		// given a vector of molecules, this will print out a PDB file
	};

}	// namespace md_files

#endif
