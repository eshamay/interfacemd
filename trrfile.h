#ifndef TRRFILE_H_
#define TRRFILE_H_

#define CPLUSPLUS 1

#include "gmx_coordinate_file.h"
#include "xdrfile/xdrfile_trr.h"

namespace gromacs {

	class TRRFile : public GMXCoordinateFile {

		public:

			TRRFile (const std::string trrpath) : GMXCoordinateFile(trrpath) { 
				// find out how many atoms are in the system
				read_trr_natoms(const_cast<char *>(this->_path.c_str()), &_natoms);

				_x = new rvec[_natoms];
				_v = new rvec[_natoms];
				_f = new rvec[_natoms];

				_coords.resize(_natoms, VecR());
				_vels.resize(_natoms, VecR());
				_forces.resize(_natoms, VecR());

				this->LoadNext();

				return; 
			}

			~TRRFile () {
				delete[] this->_x;
				delete[] this->_v;
				delete[] this->_f;
			}

			//! load the next frame of the xdr file, and parse the appropriate arrays
			void LoadNext () {

				read_trr(_file, _natoms, &_frame, &_time, &_lambda, _box, _x, _v, _f);

				// set the positions, velocities, and forces of each of the atoms
				VecR_it_non_const coord_it = _coords.begin();
				VecR_it_non_const vel_it = _vels.begin();
				VecR_it_non_const force_it = _forces.begin();

				for (int i = 0; i < _natoms; i++) {
					coord_it->Set(_x[i][0]*10.0, _x[i][1]*10.0, _x[i][2]*10.0);	// convert from nm to angstroms
					//vel_it->Set  (_v[i][0], _v[i][1], _v[i][2]);
					//force_it->Set(_f[i][0], _f[i][1], _f[i][2]);

					++coord_it; ++vel_it; ++force_it;
				}

				return;
			}

	};	// TRR File

} // namespace gromacs

#endif
