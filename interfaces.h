// Given the name of a particular residue, this function will locate the Gibb's dividing surfaces of a slab for that species. This is defined as locating the 50% point of the number density, and then returning the two locations of the interfaces.
vector<double> FindInterfaces (string atomName, string residue) {

if (_master) {
	// the technique here will be to first construct the number-density function for the residue given. This is something that the master node can do on its own.
	
	double MIN		-100.0		// minimum of the histogram
	double MAX		200.0		// maximum of the histogram
	double	DZ		0.10		// bin size

	int numbins = (int)(MAX-MIN)/DZ;
	vector<int> histo (numbins, 0);

	for (int atom = 0; atom < _numAtoms; atom++) {
		
		// find the atom/residue
		if (!((*sys)[atom]->Name() == atomName && (*sys)[atom]->Residue() == residue)) continue;	

		int bin = (int) (((*sys)[atom]->Y() - MIN) / DZ);

		histo[bin]++;
	}

	// now we have the histogram for the timestep. The next trick is to find the maxima. There are a couple scenarios:
	// 		1) The system lies within the boundaries of the box and so there are 2 interfaces
	// 		2) The system is split over the boundaries of the box, and so there will be 4 interfaces (virtually)
	// Finding the maximum, and then including all values within a certain statistical width (let's just say +/- 10%) should give us a good average value for the bulk. At that point we should be able to then find the 50% marks of the bulk values.
	
	double maxDensity = *max_element(histo.begin(), histo.end());

	double avgBulkDensity = 0.0;
	int num = 0;

	for (int i = 0; i < numbins; i++) {
		if (histo[i] > maxDensity * 0.75) {		// count only spots where the density is within 25% of the max
			num++;
			avgBulkDensity += histo[i];
		}
	}	

	avgBulkDensity /= num;

	vector<double> interfaces;
	for (int i = 0; i < numbins; i++) {
		if (histo[i] > avgBulkDensity * 0.49 && histo[i] < avgBulkDensity * 0.51)
			interfaces.push_back (DZ * i + MIN);
	}

}
		
return interfaces;
}
