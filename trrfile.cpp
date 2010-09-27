#include "trrfile.h"

	namespace gromacs {
		TRRFile::TRRFile (const std::string trrpath)
			: 
				_path(trrpath),
				_file(xdrfile_open(trrpath.c_str(), "r")) 
		{ 
			// find out how many atoms are in the system
			read_trr_natoms(const_cast<char *>(_path.c_str()), &_natoms);

			_x = new rvec[_natoms];
			_v = new rvec[_natoms];
			_f = new rvec[_natoms];

			_coords.resize(_natoms, VecR());
			_vels.resize(_natoms, VecR());
			_forces.resize(_natoms, VecR());

			this->LoadNext();

			return; 
		}

		TRRFile::~TRRFile () {
				delete[] _x;
				delete[] _v;
				delete[] _f;
				xdrfile_close(_file);
				return;
			}

		void TRRFile::LoadFirst () {
			xdrfile_close(_file);
			_file = xdrfile_open(_path.c_str(), "r");
			this->LoadNext();

			return;
		}

		// load the next frame of the trr file - coordinates, velocities, forces
		void TRRFile::LoadNext () {
			// load the frame using the xdrfile framework
			read_trr(
					_file, _natoms, &_frame, &_time, &_lambda, _box, _x, _v, _f);

			// set the positions, velocities, and forces of each of the atoms
			VecR_it_non_const coord_it = _coords.begin();
			VecR_it_non_const vel_it = _vels.begin();
			VecR_it_non_const force_it = _forces.begin();

			for (int i = 0; i < _natoms; i++) {
				coord_it->Set(_x[i][0], _x[i][1], _x[i][2]);
				vel_it->Set(_v[i][x], _v[i][y], _v[i][z]);
				force_it->Set(_f[i][x], _f[i][y], _f[i][z]);

				++coord_it; ++vel_it; ++force_it;
			}

			return;
		}

		void TRRFile::PrintInfo () const {
			printf ("Number of atoms = %d\n", _natoms);
			printf ("Frame #%d\n", _frame);
			printf ("Time = % 8.3f ps\n", _time);
			printf ("Lambda = % 8.3f\n", _lambda);

			printf ("\nBox:\n");
			this->PrintBox();

			printf ("\nCoordinates: \n");
			this->Print(_coords);
			printf ("\n");

			printf ("\nVelocities: \n");
			this->Print(_vels);
			printf ("\n");

			printf ("\nForces: \n");
			this->Print(_forces);
			printf ("\n");

			return;
		}

		void TRRFile::PrintBox () const {

			for (int i = 0; i < DIM; i++) {
				for (int j = 0; j < DIM; j++) {
					printf ("% 13.4f", _box[i][j]);
				}
				printf ("\n");
			}

			return;
		}

		void TRRFile::Print (const VecR_vec& vec) const {
			for (VecR_it it = vec.begin(); it != vec.end(); it++) {
				it->Print();
			}
			printf ("\n");
			return;
		}

	} // namespace gromacs

/*
	 int main () {
	 TRRFile trr ("sw_md.trr");
	 trr.PrintInfo();
	 trr.LoadNext();
	 trr.PrintInfo();
	 return 0;
	 }

*/
