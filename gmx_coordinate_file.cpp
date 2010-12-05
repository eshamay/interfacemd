#include "gmx_coordinate_file.h"

namespace gromacs {

	GMXCoordinateFile::GMXCoordinateFile (const std::string path)
		: 
			_path(path),
			// open the trajectory file
			_file(xdrfile_open(path.c_str(), "r")) { 
			}


	GMXCoordinateFile::~GMXCoordinateFile () {
		xdrfile_close(_file);
		return;
	}


	void GMXCoordinateFile::LoadFirst () {
		xdrfile_close(_file);
		_file = xdrfile_open(_path.c_str(), "r");
		this->LoadNext();

		return;
	}

	void GMXCoordinateFile::PrintInfo () const {
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

	void GMXCoordinateFile::PrintBox () const {

		for (int i = 0; i < DIM; i++) {
			for (int j = 0; j < DIM; j++) {
				printf ("% 13.4f", _box[i][j]);
			}
			printf ("\n");
		}

		return;
	}

	void GMXCoordinateFile::Print (const VecR_vec& vec) const {
		for (VecR_it it = vec.begin(); it != vec.end(); it++) {
			it->Print();
		}
		printf ("\n");
		return;
	}

}
