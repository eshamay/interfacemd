#pragma once
#ifndef TRRFILE_H_
#define TRRFILE_H_

#define CPLUSPLUS 1

#include "vecr.h"
#include "xdrfile/xdrfile.h"
#include "xdrfile/xdrfile_trr.h"
#include "mdfiles.h"
#include "mdsystem.h"
#include <iostream>

namespace gromacs {

	class TRRFile {

		public:

			TRRFile (const std::string trrpath);
			~TRRFile ();

			int size () const { return _natoms; }

			void LoadFirst();
			void LoadNext();

			void PrintInfo () const;

			VecR_it begin_coords () const { return _coords.begin(); }
			VecR_it end_coords () const { return _coords.end(); }

			VecR_it begin_vels () const { return _vels.begin(); }
			VecR_it end_vels () const { return _vels.end(); }

			VecR_it begin_forces () const { return _forces.begin(); }
			VecR_it end_forces () const { return _forces.end(); }

		private:

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
