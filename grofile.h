#ifndef GROFILE_H_
#define GROFILE_H_

#ifndef CPLUSPLUS
#define CPLUSPLUS 1 
#endif

#include "atom.h"
#include "molecule.h"
#include "mdfiles.h"
#include "mdsystem.h"
#include <string>
#include <vector>
#include <iostream>
#include <boost/algorithm/string.hpp>

namespace gromacs {

	class GROFile : public md_files::CoordinateFile, public MDSystem {

		public:
			GROFile (const std::string gropath);
			~GROFile ();

			void LoadFirst () { }
			void LoadNext () { }

			void Print () const;

		private:
			std::string _title;
			void _ParseSystem ();
			void _ParseMolecules () { }
	};

} // namespace gromacs

#endif
