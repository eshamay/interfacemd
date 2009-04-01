#include "mpisfg-frontend.h"
/* 
 * This script is a first attempt at dividing and conquering the analysis tasks that can take place on an amber system. For example - calculating the sum frequency spectrum of a water system requires the calculation of each individual water molecule's hyperpolarizability spectrum. Thus, we use here MPI in order to divide out the workloads.
 * Some points of note:
 * 		o Only one process is given the job of system I/O and reading/writing files. This greatly reduces the system's workload
 * 		o The entire system->molecule->atom->position/force data hierarchy utilizes the non-MPI based objects by creating a method to break down the object data, pass it around to each process, and then reconstruct it.
 * 		o This was a lot of fun. Let's do it again sometime!
 *
 */


#include <iostream>
#include <complex>

#define SYSTEM_AMBER
#include "mpisys.h"

//#include "ambermpisys.h"
#include "watersfg.h"

using std::complex;



// routine to setup the actual processing of spectral data
vector< complex<double> > ProcessData (MPIMolSystem& sys, SFGWaterAnalyzer& a);

// a way to collect all the data onto the master node
void CollectData (MPIMolSystem& sys, vector< complex<double> >& stepOutput, vector< complex<double> >& totalOutput);

/*
 * *************************************************************************************************
 */
int main (int argc, char **argv) {
	
// First things first - each node will need its own MPI system
	MPIMolSystem sys (&argc, &argv, argv[1], argv[2], argv[3]);

// run through the timesteps and calculate the hyperpolarizability spectrum for each molecule and sum them
// 	1) setup the analysis and output data structures
// 	2) master node distributes the data
// 	3) each node processes the data for its piece
// 	4) master node collects the data for the timestep
// 	5) master node advances to next timestep
	SFGWaterAnalyzer a ("SSP");	// each process gets its own sfg water analyzer
	vector< complex<double> > stepOutput;
	vector< complex<double> > totalOutput;
	totalOutput.clear();
	
	int timesteps = atoi(argv[4]);
	for (int step = 0; step < timesteps; step++) {
		
		// each node processes the molecules assigned to it
		stepOutput = ProcessData (sys, a);

		// then all the nodes congeal the data together for averaging/summing data
		if (totalOutput.empty()) totalOutput.resize(stepOutput.size(), complex<double> (0.0, 0.0));
		CollectData(sys, stepOutput, totalOutput);
		
		//printf ("%d)  %f\t%f\t%f\n", mpisys.id, waterdata.atoms[0][x], waterdata.atoms[0].Y(), waterdata.atoms[0].Z());

		// send position and force data to each node - master node broadcast, also advance to the next timestep
		sys.LoadNext();
	}

// Data output stage
	if (sys.Master()) {
		for (int i = 0; i < totalOutput.size(); i++) {
			printf ("%f\t% 13.7f\n", START_FREQ + i*FREQ_STEP, abs(totalOutput[i]));
		}
	}

return 0;
}

/*
 * *************************************************************************************************
 */

vector< complex<double> > ProcessData (MPIMolSystem& sys, SFGWaterAnalyzer& a) {

		vector< complex<double> > output;
		output.clear();

		// divide up the molecules that each process checks
		for (int mol = sys.BlockLow(sys.NumMols()); mol < sys.BlockHigh(sys.NumMols()); mol++) 
		{
			// find a water molecule
			if (sys.Molecules(mol)->Name() != "h2o") continue;
			//if (sys.Molecules(mol)->Atoms(0)->Y() < 26.0) continue;

			Water * water = static_cast<Water *>(sys.Molecules(mol));
			//(*water)["O"]->Force().Print();
			vector< complex<double> > temp;
			// process its spectrum
			temp = a.BetaSpectrum (*water);

			if (output.empty()) {	// for the first round we have to set the output size...
				output.resize(temp.size(), complex<double> (0.0, 0.0));
			}
	
			// and then add each water's spectrum onto the output
			for (int i = 0; i < temp.size(); i++) {
				output[i] += temp[i];
			}
		}

return (output);
}

void CollectData (MPIMolSystem& sys, vector< complex<double> >& stepOutput, vector< complex<double> >& totalOutput) {

	// this is where the output of each process is pulled together on the master node. However, the complex data has to be broken down into its real and imaginary pieces first, and then sent, before sticking it all back together on the master node
	// So first let's get each node to break down the data for shipping it out
	double real[stepOutput.size()];
	double imag[stepOutput.size()];

	double realdata[stepOutput.size()];
	double imagdata[stepOutput.size()];

	// now get all the data into a nice spot in memory to be reduced to the front node
	for (int i = 0; i < stepOutput.size(); i++) {
		real[i] = stepOutput[i].real();
		realdata[i] = 0.0;
		imag[i] = stepOutput[i].imag();
		imagdata[i] = 0.0;
	}

	// let's get all the processes together now
	//MPI_Barrier (sys.WorldComm());

	// grab the sum of all the data computed until now
	// here's the real data
	MPI_Reduce (real, realdata, stepOutput.size(), MPI::DOUBLE, MPI::SUM, 0, sys.WorldComm());

	// and here's the imaginary data
	MPI_Reduce (imag, imagdata, stepOutput.size(), MPI::DOUBLE, MPI::SUM, 0, sys.WorldComm());
	
	if (sys.Master()) {
		// now run through the data and add it on to the master node's copy
		for (int i = 0; i < stepOutput.size(); i++) {
			complex<double> data (realdata[i], imagdata[i]);
			totalOutput[i] += data;
		}
	}
	
	/*
	if (!sys.Master()) {
	
		for (int i = 0; i < output.size(); i++) {
			real[i] = output[i].real();
			imag[i] = output[i].imag();
		}
	
		// now each node waits for the master node to ask for the data before shipping it out
		MPI_Recv (&msg, 1, MPI_INT, 0, handshaketag, sys.WorldComm(), sys.Stat());
			printf("%d listening for ping\n", sys.ID());

		// Once the master has asked for the data, send it out
		MPI_Send (real, output.size(), MPI_DOUBLE, 0, xfertag, sys.WorldComm());
			printf("%d sends real data\n", sys.ID());
		// then wait for the master to say everything is cool
		MPI_Recv (&msg, 1, MPI_INT, 0, handshaketag, sys.WorldComm(), sys.Stat());
			printf("%d listening for ping\n", sys.ID());
		// then send the imaginary data
		MPI_Send (imag, output.size(), MPI_DOUBLE, 0, xfertag, sys.WorldComm());
			printf("%d sends imag data\n", sys.ID());
	}

	if (sys.Master()) {
		
		// the master node goes through each process and sends out a handshake call to ask for data
		for (int proc = 1; proc < sys.Procs(); proc++) {
			MPI_Send (&msg, 1, MPI_INT, proc, handshaketag, sys.WorldComm());
				printf("master pings %d\n", proc);

			// once the processes have shook hands, the master receives the data
			MPI_Recv (real, output.size(), MPI_DOUBLE, proc, xfertag, sys.WorldComm(), sys.Stat());
				printf("master receiving real data from %d\n", proc);
			// and then when it's been processed, let's ask again for the imaginary data
			MPI_Send (&msg, 1, MPI_INT, proc, handshaketag, sys.WorldComm());
				printf("master pings %d\n", proc);
			// and then get that! that's that!
			MPI_Recv (imag, output.size(), MPI_DOUBLE, proc, xfertag, sys.WorldComm(), sys.Stat());
				printf("master receiving real data from %d\n", proc);

			// now run through the data and add it on to the master node's copy
			for (int i = 0; i < output.size(); i++) {
				output[i] += complex<double> (real[i], imag[i]);
			}
		}
	}
	*/
return;
}
