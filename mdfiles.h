#ifndef MD_FILES_H_
#define MD_FILES_H_

#include "atom.h"
#include "molecule.h"
#include <iostream>

namespace md_files {

	class CoordinateFile {

		protected:
			FILE *_file;				// the XYZ file listing all the atom coordinates
			std::string _path;

			char _line[1000];

		public:

			typedef Eigen::Map<VecR>	coord_t;
			typedef std::vector<coord_t> coord_set_t;
			typedef coord_set_t::const_iterator coord_it;

			CoordinateFile (const std::string path) 
				:
					_file ((FILE *)NULL),
					_path(path) {

						_file = fopen (path.c_str(), "r");
						if (_file == (FILE *)NULL)
						{
							printf ("Couldn't load the Coordinate file %s\n", path.c_str());
							exit(1);
						}
					}


			virtual ~CoordinateFile () {
				fclose (_file);
			}

			char * ReadLine () { 
				fgets (_line, 1000, _file);
				return _line;
			}

			char * Line () {
				return _line;
			}

			virtual void LoadNext () = 0;

	};	// Coordinate file


}	// namespace md_files

#endif
