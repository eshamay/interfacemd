#include "mpitmp.h"
#include <iostream>
#include <complex>
#include "mpisys.h"
#include "dipoleparm.h"

using std::complex;


int main (int argc, char **argv) {

	string xyzpath = "/home/users/eshamay/MD/mdcrd/waterslab/traj.xyz";
	VecR syssize (26.0, 52.0, 26.0);
	string wannierpath = "";
	MPIMolSystem sys (&argc, &argv, xyzpath, syssize, wannierpath);
	WaterDipoleParms dips ("/home/users/eshamay/work/Water/dipole-TCF-spectra/water.dipole.parms.dat");   // file containing all the dipole info

	double center;
	Water * wat;

	for (int step = 0; step < 50000; step++) {

		double totalDip[3] = {0.0, 0.0, 0.0};

		if (sys.Master()) {

			vector<double> position;
			position.clear();
	
			// let's first find the extent of molecules on the long axis to determine where the system-center is located
	    	for (int atom = 0; atom < sys.NumAtoms(); atom++) {
				position.push_back (sys.Atoms(atom)->Y());
			}
			sort(position.begin(), position.end());
	
			double min = position[0];
			double max = position[position.size()-1];
			center = (max - min)/2.0;
		}
	
		MPI_Bcast (&center, 1, MPI_DOUBLE, 0, sys.WorldComm());

		// now each node processes through a chunk of the molecules

		for (int mol = sys.BlockLow(sys.NumMols()); mol < sys.BlockHigh(sys.NumMols()); mol++) {
			if (sys.Molecules(mol)->Name() != "h2o") continue;

			wat = static_cast<Water *>(sys.Molecules(mol));
			wat->FindMolecularAxes();

			VecR dip;
			dip = dips.Dipole (wat->OH1()->Magnitude(), wat->OH2()->Magnitude(), acos((*wat->OH1()) < (*wat->OH2()))*180/M_PI);
			wat->RotateToLab(dip);
			if (wat->Atoms()[0]->Y() < center) {
				dip.Y(dip[y] * -1.0);
			}
			// update the total dipole
			totalDip[0] += dip[x];
			totalDip[1] += dip[y];
			totalDip[2] += dip[z];
		}

		double dipX, dipY, dipZ;
		// then we gather all the data together
		MPI_Allreduce (&totalDip[0], &dipX, 1, MPI_DOUBLE, MPI_SUM, sys.WorldComm());
		MPI_Allreduce (&totalDip[1], &dipY, 1, MPI_DOUBLE, MPI_SUM, sys.WorldComm());
		MPI_Allreduce (&totalDip[2], &dipZ, 1, MPI_DOUBLE, MPI_SUM, sys.WorldComm());

		if (sys.Master()) {
			printf ("% 8d\t% 12.8f\t% 12.8f\t% 12.8f\n", step, dipX, dipY, dipZ); 
			fflush(stdout);
			sys.LoadNext();
		}
	}

return 0;
}
