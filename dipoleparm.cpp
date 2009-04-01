#include "dipoleparm.h"

WaterDipoleParms::WaterDipoleParms (string parmpath="dipoleparm.dat") {
	
	// firstly, we need to find out some information about the file's contents. We assume here that the file has a single line of header information, and the rest is data layed out as follows:
	// 		R1	R2	Theta	X	Y	Z	Norm
	//
	// 	Where -
	// 		R1		= Length of OH bond 1
	// 		R2		= Length of OH bond 2
	// 		Theta	= Angle of the water molecule
	// 		X, Y, Z	= vector of the dipole moment
	// 		Norm	= Norm (magnitude) of the dipole moment vector
	// 
	// So let's find out about the file and load it's contents into some data structures!
	
	_file = (FILE *)NULL;
	_file = fopen(parmpath.c_str(), "r");
	if (_file == (FILE *)NULL) {
		printf ("Sorry, the parmfile %s doesn't exist\n", parmpath.c_str());
		exit (1);
	}

	// Now we have the task of finding out how many data points we have for each parameter (R1, R2, and theta). We'll use some shell commands to dig out the data
	
	for (int i = 0; i < 3; i++) {
		// here's the commands we'll use to get some info
		char col [10];
		sprintf (col, "%d", i+1);
		string col_str=col;

		// this should get us the number of bins for each parameter
		string bin_com = "awk '{print $" + col_str + "}' " + parmpath + " | sort | uniq | wc -l";
		FILE * fp = popen(bin_com.c_str(), "r");
		fscanf (fp, "%d", &_num_bins[i]);
		_num_bins[i]--;
		pclose(fp);
		
		// to get the max values
		string max_com = "awk '{print $" + col_str + "}' " + parmpath + " | sort -g | uniq | tail -n 1";
		fp = popen(max_com.c_str(), "r");
		fscanf (fp, "%lf", &_max[i]);
		pclose(fp);

		// to get the min values
		string min_com = "awk '{print $" + col_str + "}' " + parmpath + " | sort -g | uniq | head -n 2 | tail -n 1";
		fp = popen(min_com.c_str(), "r");
		fscanf (fp, "%lf", &_min[i]);
		pclose(fp);

		// let's also set the bin size for each parameter
		_dr[i] = (_max[i] - _min[i]) / (_num_bins[i] - 1);
	}

	// here's where we set up all the pointer arrays to the final data values of the dipole vector coords and the norm
	_data = (double ****) malloc (_num_bins[0] * sizeof(double ***));

	for (int r1 = 0; r1 < _num_bins[0]; r1++) {
		_data[r1] = (double ***) malloc (_num_bins[1] * sizeof(double **));

		for (int r2 = 0; r2 < _num_bins[1]; r2++) {
			_data[r1][r2] = (double **) malloc (_num_bins[2] * sizeof(double *));

			for (int theta = 0; theta < _num_bins[2]; theta++) {
				_data[r1][r2][theta] = (double *) malloc (4 * sizeof(double));
			}
		}
	}	

		
	// the final step of setup is to read through the file and pull each bit of data into the _data holder
	rewind(_file);
	fscanf (_file, " %*s %*s %*s %*s %*s %*s %*s");		// skip the header

	double r1, r2, theta, x, y, z, norm;
	double * data;

	while (!feof(_file)) {
		fscanf (_file, " %lf %lf %lf %lf %lf %lf %lf", &r1, &r2, &theta, &x, &y, &z, &norm);		// grab the data from the file
		//printf ("%f\t%f\t%f\n", x, y, z);
		data = _Data (r1, r2, theta);
		data[0] = x;
		data[1] = y;
		data[2] = z;
		data[3] = norm;
	}
	
return;
}

// this will run through and pull out the final data location for a given set of parameters
double * WaterDipoleParms::_Data (double r1, double r2, double theta) {

	// first some safety checks to make sure that the parameters are within the limits
	if (r1 < _min[0] || r1 > _max[0] + _dr[0] \
		|| r2 < _min[1] || r2 > _max[1] + _dr[1] \
		|| theta < _min[2] || theta > _max[2]) {

		printf("dipoleparm.h) _Data() :: parameters are out of limits, check r1, r2, and theta values.\n%f\t%f\t%f\n", r1, r2, theta);
		exit(1);
	}

	// now let's determine the bin number for each of the parameters
	int r1_bin = (int)((r1+.000001-_min[0])/_dr[0]);
	int r2_bin = (int)((r2+.000001-_min[1])/_dr[1]);
	int theta_bin = (int)((theta+.000001-_min[2])/_dr[2]);

	//printf ("(%f\t%f\t%f)\t%d %d %d\n", r1, r2, theta, r1_bin, r2_bin, theta_bin);
	//printf ("%d %d %d\n", r1_bin, r2_bin, theta_bin);
return _data[r1_bin][r2_bin][theta_bin];
}
