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

	};	// Coordinate file


}	// namespace md_files

#endif
