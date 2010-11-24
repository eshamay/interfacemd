#ifndef GMX_COORDINATE_FILE_H_
#define GMX_COORDINATE_FILE_H_

#ifndef CPLUSCPLUS
#define CPLUSPLUS 1
#endif

#include "vecr.h"
#include "xdrfile/xdrfile.h"
#include <iostream>

namespace gromacs {

	class GMXCoordinateFile {

		public:

			GMXCoordinateFile (const std::string path);
			virtual ~GMXCoordinateFile ();

			int size () const { return _natoms; }

			virtual void LoadFirst();
			//! load the next frame of the xdr file, and parse the appropriate arrays
			virtual void LoadNext() = 0;

			//! Prints a whole lot of information about the xdr/coordinate file (e.g. box size, number of atoms, etc)
			virtual void PrintInfo () const;

			VecR_it begin_coords () const { return _coords.begin(); }
			VecR_it end_coords () const { return _coords.end(); }

			VecR_it begin_vels () const { return _vels.begin(); }
			VecR_it end_vels () const { return _vels.end(); }

			VecR_it begin_forces () const { return _forces.begin(); }
			VecR_it end_forces () const { return _forces.end(); }

		protected:

			std::string _path;
			XDRFILE * _file;

			int	_natoms;

			int _frame;
			float _time;
			float _lambda;
			matrix _box;
			rvec *_x, *_v, *_f;	// coordinates, velocities, and forces from each timeframe
			VecR_vec _coords, _vels, _forces;

			void PrintBox () const;
			void Print (const VecR_vec& vec) const;
	};



} // namespace gromacs

#endif
