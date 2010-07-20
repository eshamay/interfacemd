/* These routines were originally written for analysis of nitric acid, with the .xyz file format specifically used. */
#include "watersystem.h"
#include <iostream>
#include <string>
#include <math.h>
#include <stdio.h>


#define NO_BOND	2.0

typedef vector<double> trajectory;

class WaterAnalyzer {

	WaterSystem&			_sys;		// the water system that will be used for calculations

public:

	WaterAnalyzer (WaterSystem& system) : _sys(system) {}		// Constructor that grabs the system and loads it up

	WaterSystem&	System() { return _sys; }

	vector< vector<int> > ZCoordHistogram ( double DZ );	// Create a 2D histogram of Z-axis position and h-bond coordination for waters
	trajectory DipoleCorrelation ();						// Create a TCF of the system polarization
	trajectory SysDipoleCorrelation ();						// Create a TCF of the system polarization
	vector<double> BondTrajectory (string atomName1, string atomName2, double bondlengthmax);
	vector<double> RDF (string atomName1, string atomName2);	// Calculate an RDF from the system
	vector< vector<trajectory> > CoordTrajectories ();			// Finds all OH-bond trajectories for a all coordinations of water
	trajectory NO3Rotation ();
};


/*
Here we create a 2D histogram that gives intensities for the z-axis position of a water molecule, and the coordination of that
same molecule. The vector created has 2 indices: [z-position][coordination]. The value stored is the intensity (particle density) over
the entire timespan of the system.
*/
vector< vector<int> > WaterAnalyzer::ZCoordHistogram ( double DZ ) {

	// the z-axis is now cut up into pieces sized DZ wide (i.e. each bin covers an amount DZ on the axis)
	// Then each coordination is accounted for at each z-position by using a vector to hold an element for each coord type
	// here it's 12 = number of different coord types (as defined in the watersystem.h enum)
	vector< vector<int> > histogram ((int)(_sys.Size()[z] / DZ), vector<int>(12, 0));
	int zbin;
	coordination coord;

	// prep the system to the first timestep
	_sys.LoadFirst();

	// run through each time step and bin the values
	while (_sys.Current() != _sys.Last()) {
		_sys.FindCoordinations();

		// note: water position is based on the oxygen's z-position; not the center of mass
		RUN (_sys.Oxygens()) {
			zbin = (int)(_sys[_sys.Oxygens()[i]][z] / DZ);
			coord = _sys.Coordinations()[i];

			histogram[zbin][coord]++;		// update the correct bin
		}

		// and lastly move to the next timestep
		_sys.LoadNext();
	}

return (histogram);
}

/*
A routine for calculating the dipole-dipole autocorrelation function for a given trajectory.
As per Meuwly (2003), a slice of the trajectory, N timesteps, is used. N is chosen such that 2^N comprises 1/3 to 1/2 of the total timesteps.
A vector with N elements of the trajectory is then produced beginning with the first timestep and the correlation function is calculated.
Then the trajectory is shifted to i+1 and so on until all the trajectory has been used. At this point the correlations are averaged and shipped off as a result.
*/
trajectory WaterAnalyzer::DipoleCorrelation () {

	// First, instead of having to recalculate the system dipoles EACH time we load a new timestep... and since the dipoles are the only
	// quantity we're interested in, let's make a 2D vector that will hold the dipole of each molecule at each timestep.
	// This means running through the entire trajectory once in order to grab all the dipole vectors of each molecule at each timestep.
	_sys.LoadFirst();
	vector< vector<VecR> > dipoles;				// this will be our dipole [timestep][molecule] vector
	while (_sys.Current() != _sys.Last()) {
		_sys.FindDipoles();
		dipoles.push_back(_sys.Dipoles());
		_sys.LoadNext();
	}

	// okay, the bulky time-consuming stuff is mostly done. Let's get on with the TCF calculations
	bool first = true;		// are we in the first timestep of each running-chunk of the trajectory?

	int N  = (int)((_sys.Last() - _sys.First()) / 4);		// the number of steps in the moving correlation queue for processing

	int start = 0, 							// the starting timestep of each TCF chunk
		step,								// current timestep being processed
		last = _sys.Last() - _sys.First();	// the last timestep available

	// here we advance our time chunk over the entire trajectory until our chunk hits the ending timestep.
	vector<trajectory> chunkTCFstack;
	while (start < last - N) {

		// this will hold all the TCF's of the individual waters for a given trajectory chunk
		// TCFstack [water][timestep]
		vector<trajectory> TCFstack (dipoles[0].size(), trajectory());

		// then run through each timestep from (start) to (start + N) and load the correlation for each molecule
		// i.e. go through each timestep in the trajectory chunk
		for (step = start; step < start+N; step++) {

			// for each timestep we'll calculate the correlation (dot-product) of each water's dipole with that of time t=0
			RUN (dipoles[step]) {

				if (step == start) {		// for the first step in each trajectory chunk, just record the dipole for time = 0
					TCFstack[i].clear();	// just in case there's some residue in our vector
					VecR vZeroDipole = dipoles[start][i];
					//TCFstack[i].push_back( (vZeroDipole * vZeroDipole) / vZeroDipole.Magnitude() / vZeroDipole.Magnitude() ); // should=1.0
					TCFstack[i].push_back( vZeroDipole[z] * vZeroDipole[z] );
				}

				else {
					// for every other step, record the correlation with the t=0 step
					// note - normalizing to the t(0) dipole magnitude stored in zeroT
					VecR vZeroDipole = dipoles[start][i];
					//TCFstack[i].push_back( (dipoles[step][i] * vZeroDipole) / vZeroDipole.Magnitude() / vZeroDipole.Magnitude() );
					TCFstack[i].push_back( (dipoles[step][i][z] * vZeroDipole[z]) );
				}
			}
		}

		/*
		At this point we have a 2D vector (TCFstack) that represents our trajectory chunk (indexed as [molecule][timestep]). Now we should try one of two things depending on our mood - we can average all the individual molecular TCFs into a single TCF that represents our TCF chunk, or we can fourier-transform each TCF and then average the spectra... For now, I'll write the code to average the TCFs, and then fourier transform the output somewhere outside, as I believe that ultimately the two methods are equivalent (but I'm 100% sure on that point!) :)
		so here goes - average all the TCF into a single chunk-TCF
		*/
		trajectory chunkTCF;
		RUN (TCFstack[0]) {	// i = timestep #
			double total = 0.0;
			RUN2 (TCFstack) {		// j = molecule
				 // for each timestep let's add up all the TCF values for each molecule
				 total += TCFstack[j][i];
			}
			total /= TCFstack.size();	// average the value over the # of molecules
			chunkTCF.push_back (total);	// and add in the value for the current timestep
		}

		/*
		Now we have a 1D trajectory that represents the averaged dipole TCF of all the molecules over that chunk of the total trajectory. SO... let's store it in some vector so that we can later average that with all the other trajectory-chunk TCFs for a final, very smooth, TCF
		*/
		chunkTCFstack.push_back (chunkTCF);

		// and then we have to advance in the total trajectory to process a new chunk
		start++;
	}

	/*
	Just like before, we have a 2D stack of TCF's, so let's average all of those like before.
	*/
	trajectory output;

	RUN (chunkTCFstack) {
		double total = 0.0;
		RUN2 (chunkTCFstack[i]) {
			total += chunkTCFstack[j][i];
		}
		total /= chunkTCFstack.size();
		output.push_back (total);
	}
return output;
}


/*	This is older code - deprecated 04 Jan 2007
		// Now we've got a 2D vector (TCF) that holds a lot of correlation information. Each element of TCF is a new timestep, and holds
		// the correlation of each molecule at that timestep. Thus the indices are TCF[timestep][molecule], and the value is what was
		// computed above (normalized correlation to the t(0) dipole of the molecule).
		// Let's compress all the molecular correlations into a
		// single value as the average molecular correlation for each timestep. The resulting 1D vector will be our "queue" that
		// represents a single correlationg function of length N pulled from somewhere in the trajectory.
		// In conclusion - compress the 2D vector of molecular correlations to a 1D of average TCF
		currentQueue.clear();
		RUN (TCF) {
			currentQueue.push_back(0.0);
			// here the molecular correlations are summed
			RUN2 (stepTCF) {
				currentQueue[i] += TCF[i][j];
			}
			// and now we average by dividing by the number of molecules probed
			currentQueue[i] = currentQueue[i] / dipoles[step].size();
		}

		// now that the single-queue TCF has been calculated for each molecule, averaged, and polished, let's shove it onto a stack of queues
		queues.push_back(currentQueue);

		//printf (" and ended on step %d\n", step);
		start++;		// And lastly, update the starting timestep so that the next queue can be calculated as we chug our way
						// through the entire trajectory
	}

	// now we've got a vector full of queues! Each queue is one correlation function sampled from somewhere else in the trajectory.
	// So now, as before, let's compress all those TCFs into a single averaged TCF that we can call an output and be done with it.
	trajectory output;
	output.clear();
	RUN (currentQueue) {
		output.push_back(0.0);
		RUN2 (queues) {
			output[i] += queues[i][j];
		}
		output[i] = output[i] / queues.size();
	}

return (output);
}
*/

trajectory WaterAnalyzer::SysDipoleCorrelation () {

	// the output time-correlation function vector
	trajectory TCF;		// A TCF pulled from some part of the total trajectory. The length of this is determined by N (below)
	vector<trajectory> TCFStack;		// all the TCFs from over the trajectory
	VecR vZero;				// The t=0 total system dipole
	VecR vStep;				// the system dipole at time t

	bool first = true;		// are we in the first timestep of each running-chunk of the trajectory?

	int N  = (int)((double)(_sys.Last() - _sys.First()) / 2.5);		// the number of steps in the moving correlation queue for processing
	//cout << N << " steps per correlation\n";

	// First, instead of having to recalculate the system dipoles EACH time we load a new timestep... and since the dipoles are the only
	// quantity we're interested in, let's make a 2D vector that will hold the dipole of each molecule at each timestep.
	// This means running through the entire trajectory once in order to grab all the dipole vectors of each molecule at each timestep.
	_sys.LoadFirst();
	vector< vector<VecR> > dipoles;				// this will be our dipole [timestep][molecule] vector
	while (_sys.Current() != _sys.Last()) {
		_sys.FindDipoles();
		dipoles.push_back(_sys.Dipoles());
		_sys.LoadNext();
	}

	// okay, the bulky time-consuming stuff is mostly done. Let's get on with the TCF calculations

	int start = 0, 							// the starting timestep of each TCF chunk
		step,								// current timestep being processed
		last = _sys.Last() - _sys.First();	// the last timestep available

	// here we advance our TCF chunk over the trajectory until our chunk hits the ending timestep.
	while (start < last - N) {
		first = true;		// reset the first-step marker to grab the t(0) molecular dipoles (loaded into the zeroT vector)
		TCF.clear();

		// run through each timestep from (start) to (start + N) and load the total system dipole correlation
		for (step = start; step < start+N; step++) {

			vStep.Zero();
			// sum all the dipoles in the system for the total system polarization
			RUN (dipoles[step]) {
				vStep += dipoles[step][i];
			}

			// if this is the start of a new TCF chunk, then set the t=0 dipole for auto-correlation calcs later down
			if (first) {
				vZero = vStep;
				first = false;
			}

			// note - normalizing to the t(0) dipole magnitude stored in vZero
			TCF.push_back ( (vStep * vZero) / vZero.Magnitude() / vZero.Magnitude() );

		}

		// now that the TCF has been calculated for each molecule, averaged, and polished, let's shove it onto a stack of queues
		TCFStack.push_back (TCF);

		start++;		// And lastly, update the starting timestep so that the next queue can be calculated as we chug our way
						// through the entire trajectory
	}

	// Now we've got a stack of TCFs. Let's average each timestep's correlation into a single output vector
	trajectory output (TCF.size(), 0.0);

	RUN (TCFStack) {
		RUN2 (TCF) {
			output[j] += TCFStack[i][j];
		}
	}

	RUN (output) {
		output[i] /= TCFStack.size();	// and average over the number of trajectories processed
	}

return (output);
}

/* Produce a plot of the bond-length between two atom species. i.e. if "O" and "H" are chosen, then all bonds between an O and an H that are shorter than the max length (bondlengthmax) are added together for each time-step, and then output. */
vector<double> WaterAnalyzer::BondTrajectory (string atomName1, string atomName2, double bondlengthmax) {

	vector<double> results;
	double stepSum;			// a running total for the current timestep
	double distance;		// caclculated bond distance between particles
	int bonds;			// number of OH bonds processed

	_sys.LoadFirst();	// rewind the system
	while (_sys.Current() != _sys.Last()) {
		stepSum = 0.0;
		bonds = 0;

		RUN (_sys) {

			// Find the occurence of each atom1
			if (_sys[i].Name() != atomName1) continue;

			// and then each atom2 for measuring the distance
			RUN2 (_sys) {
				if (_sys[j].Name() != atomName2) continue;

				// calculate the distance, see if it's within a bond-length (or whatever the max is set to)
				distance = _sys[i] - _sys[j];
				if (distance < bondlengthmax) {
					// and then add it into the trajectory sum-total
					stepSum += distance;
					bonds++;
				}
			}
		}

		if (bonds) stepSum = stepSum / (double)bonds;
		results.push_back (stepSum);

		_sys.LoadNext();
	}

return (results);
}

/* A routine for calculating a radial-distribution function, or pair correlation function. Look up your favorite text for a definition of RDFs. I personally enjoy them because they give a really concise definition of bonding distances (equilibrium and extrema) for bond-lengths. The following routine normalizes the output so that as the distance from a particle increases, the function goes to 1.0. */
vector<double> WaterAnalyzer::RDF (string atomName1, string atomName2) {

	_sys.LoadFirst();	// rewind the system

	double maxdistance = 6.5;	// The max distance to check for the RDF (largest 'shell') in angstroms
	double DR = .01;			// the differential shell radius (decrease this for higher precision... but .01 is already pretty good)
	double volume = 4.0/3.0*M_PI*pow(maxdistance,3);		// total volume over which we're checking

	vector< vector<double> > results;		// for holding the RDF of each timestep (before averaging over the # of timesteps)

	double distance;			// distance between two particles

	// let's run over each timestep
	while (_sys.Current() != _sys.Last()) {

		int num = 0;	// number of particle-pairs processed

		vector<int> histogram ((maxdistance/DR)+1, 0);
		// find each occurence of the first particle by name
		RUN (_sys) {
			if (_sys[i].Name() != atomName1) continue;

			// and also every occurence of the 2nd particle by name		- note that we don't cover particles twice!
			RUN2 (_sys) {
				if (_sys[j].Name() != atomName2) continue;

				// now calculate the distance between the two and bin into the right place in the histogram
				distance = _sys[i] - _sys[j];

				if (distance < maxdistance) {
					num++;	// and since we have a pair, let's increase our counter
					int bin = (int)(distance/DR);	// here's the bin
					//printf ("%d % .3f- %d\t%d\n", j, distance, bin, histogram[bin]);
					histogram[bin]++;		// and this is the histogram that Jack built
				}
			}
		}

		// That wraps up the number density calculations where the interparticle pair distances were binned.
		// Now we calculate the g(r) function for each differential distance based on the number of
		// particles in each differential shell.

		/* the following normalization treatment was pulled from "Essentials of Computational Chemistry (2nd Ed.) - Theories
		   and Models", Christopher J Cramer. Wiley 2004. pp. 84-86.
		*/
		double norm;		// the normalization constant for each distance (r)
		vector<double> RDF (histogram.size(), 0.0);		// an output histogram with normalization applied to each bin

		for (int i = 1; i < histogram.size(); i++) {
			norm = volume / (4 * M_PI * pow(i*DR,2) * DR) / num ;		// here's our scaling constant
			// and then normalize g(r)
			RDF[i] = histogram[i] * norm;
		}

		results.push_back (RDF);

		_sys.LoadNext();
	}
	// lastly, let's average through all the RDFs we've collected
	vector<double> output;

	RUN (results[0]) {
		output.push_back(0.0);
		RUN2 (results) {
			output[i] += results[j][i];
		}
		output[i] = output[i] / results.size();
	}

	return (output);
}

/* This is an attempt to generate spectra for each individual coordination type of water molecules. The routine is rather obtuse, but I'll attempt to explain what's going on. First, what's expected:
	A single water molecule will remain in a given coordination for a finite amount of time before on h-bond breaks or forms, and thus changes its coordination. For this reason, if we want a spectrum for molecules with a given coordination, we cannot just grab all waters that start out in that coordination and monitor them over the trajectory; we must, therefor, pick out the parts of the trajectory for each molecule where it sits in the given coordination.
	An issue arises when we consider that a "good" spectrum requires a certain amount of timesteps, (i.e. we can't generate a meaningful spectrum from 50 timesteps) so what I've done here is a somewhat obfuscated system of building trajectory "chunks" for each water molecule over the trajectory for which they sit in the desired coordination. The chunks are only considered valid if they are a minimum number of timesteps (500?). If they are atleast this big, then the chunk is added to a pile of chunks of the trajectory for processing.
	Once all the pieces of the trajectory have been gathered for each coordination type, they are fourier-transformed into spectra, and then compiled into a single spectrum.
*/

#define MIN_TRAJ_SIZE	600
vector< vector<trajectory> > WaterAnalyzer::CoordTrajectories () {

	vector< vector<trajectory> > pile (12, vector<trajectory>());	// a pile to hold all the trajectory-chunks for each coordination

	vector< trajectory > assembly (_sys.size(), trajectory());	// where the trajectory-chunks will be built for each molecule

	/* the markers vector is sized according to the system-size because the sys.Oxygens().size() may change, because in any given frame it's possible that the # waters will change (to OH- or H3O+)
	the markers keep track of the current and previous coordination-type of each water molecule. This is how we know when the coordination changes and when we need to toss the gathered trajectory-chunk on the growin pile for processing, or trash it because it's too short.
	*/
	vector <coordination> markers (_sys.size(), UNBOUND);

	// let's set the first coordination settings/initial conditions
	_sys.LoadFirst();
	_sys.FindCoordinations();
	RUN (_sys.Oxygens()) {
		markers[_sys.Oxygens()[i]] = _sys.Coordinations()[i];
	}

	while (_sys.Current() != _sys.Last()) {		// for each timestep...

		RUN (_sys.Oxygens()) {								// for each water...
			int oxy = _sys.Oxygens()[i];					// find the water's index in the system
			coordination coord = _sys.Coordinations()[i];	// and its coordination

			// If the coordination is not the same as the previous timestep then the water has changed coordination. At this point we have to check to see if the trajectory-chunk has built to a large enough size and whether it should be processed or not. If it is big enough to be useful, then do a little magic with finishing off the trajectory chunk and putting it on the growing pile (sorted by coordination)
			if (coord != markers[oxy]) {
				// check to see that it is long enough for a decent spectrum (# of timesteps??)
				if (assembly[oxy].size() >= MIN_TRAJ_SIZE) {
					// toss it on the pile (in the correct coordination bin)
					pile[coord].push_back(assembly[oxy]);
					// and clear out the molecule's spot in the assembly line for a new trajectory before starting to build it up
					assembly[oxy].clear();
				}
			}

			// Let's now turn our attention to building up chunks of the trajectory
			double bondsum = 0.0;
			RUN2(_sys.Map()) {
				if (_sys.Map()[oxy][j] > 0.0 && oxy != j) {	// run through the waterMap connectivity matrix and find each OH bond
					bondsum += _sys.Map()[oxy][j];			// sum the two values of the OH-bond length
				}
			}
			assembly[oxy].push_back(bondsum/2.0);	// then add it to the assembly line and average between the two bonds

			// and update the coordination marker for the next frame
			markers[oxy] = coord;
		}

		// lastly, update to the next frame of the trajectory, and find all the new coordinations
		_sys.LoadNext();
		_sys.FindCoordinations();
	}

return (pile);
}

/* find the nitric acid molecule by locating the Nitrogen and 3 closest oxygens (within an NO bond-length). Then calculate the rotation wrt the z-axis. */
trajectory WaterAnalyzer::NO3Rotation () {

	trajectory output;

	// indices for the atoms in the xyz file
	int N = 0,
		O1 = 0,
		O2 = 0,
		O3 = 0;

	// find the index of each atom in the nitric acid
	RUN (_sys) {
		if (_sys[i].Name() == "N") {
			N = i;
		}
	}
	RUN (_sys) {
		if (_sys[i].Name() == "O" && _sys[N] - _sys[i] < NO_BOND) {
			if (!O1) O1 = i;
			else if (!O2) O2 = i;
			else if (!O3) O3 = i;
		}
	}

	VecR Z (0.0, 0.0, 1.0);

	while (_sys.Current() != _sys.Last()) {

		VecR NO1 = _sys[N].Position() - _sys[O2].Position();
		VecR NO2 = _sys[N].Position() - _sys[O3].Position();
		VecR Normal = NO1 % NO2;
		double value = Normal < Z;
		output.push_back(value);

		_sys.LoadNext();
	}

return output;
}
