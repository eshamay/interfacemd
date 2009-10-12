
/* A routine for calculating a radial-distribution function, or pair correlation function. Look up your favorite text for a definition of RDFs. I personally enjoy them because they give a really concise definition of bonding distances (equilibrium and extrema) for bond-lengths. The following routine normalizes the output so that as the distance from a particle increases, the function goes to 1.0. */
vector<double> RDF (string atomName1, string atomName2, vector<Atom *>& sys) {

#define maxdistance 15.0	// The max distance to check for the RDF (largest 'shell') in angstroms
#define DR	.01			// the differential shell radius (decrease this for higher precision... but .01 is already pretty good)

  double volume = 4.0/3.0*M_PI*pow(maxdistance,3);		// total volume over which we're checking
  double distance;			// distance between two particles
  int num = 0;				// number of particle-pairs processed

  vector<int> histogram ((maxdistance/DR)+1, 0);
  // find each occurence of the first particle by name
  RUN (sys) {
	if (sys[i]->Name().find(atomName1) == string::npos) continue;

	// and also every occurence of the 2nd particle by name		- note that we don't cover particles twice!
	RUN2 (sys) {
	  if (sys[j]->Name().find(atomName2) == string::npos) continue;
	  /* Don't look at the radial-distance of an atom to itself */
	  if (sys[i] == sys[j]) continue;
	  // now calculate the distance between the two and bin into the right place in the histogram
	  distance = sys[i]->MinDistance(sys[j]);

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

  RUN (histogram) {
	
	norm = volume / (4 * M_PI * pow(i*DR,2) * DR) / num ;		// here's our scaling constant
	// and then normalize g(r)
	RDF[i] = histogram[i] * norm;
  }

  // lastly, let's average through all the RDFs we've collected
  vector<double> output;

  return (histogram);
}
