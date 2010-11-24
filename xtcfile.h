#ifndef XTCFILE_H_
#define XTCFILE_H_

#define CPLUSPLUS 1

//#include "gmx_coordinate_file.h"
#include "xdrfile/xdrfile_xtc.h"

namespace gromacs {

	class XTCFile : public GMXCoordinateFile {

		public:

			XTCFile (const std::string xtcpath) : GMXCoordinateFile (xtcpath) { 

				// find out how many atoms are in the system
				read_xtc_natoms(const_cast<char *>(this->_path.c_str()), &this->_natoms);

				this->_x = new rvec[_natoms];

				this->_coords.resize(_natoms, VecR());
				this->_vels.resize(_natoms, VecR());
				this->_forces.resize(_natoms, VecR());

				this->LoadNext();

				return; 
			}

			~XTCFile () {
				delete [] this->_x;
			}

			void LoadNext () {
				read_xtc(this->_file, this->_natoms, &this->_frame, &this->_time, this->_box, this->_x, &this->_prec);
				// set only the positions because xtc files don't have anything else
				VecR_it_non_const coord_it = _coords.begin();

				for (int i = 0; i < _natoms; i++) {
					coord_it->Set(_x[i][0]*10.0, _x[i][1]*10.0, _x[i][2]*10.0);	// convert from nm to angstroms

					++coord_it;
				}
				return;
			}


		private:
			float _prec;	// precision that the xtc was stored in
	};	// XTC File

} // namespace gromacs

#endif
