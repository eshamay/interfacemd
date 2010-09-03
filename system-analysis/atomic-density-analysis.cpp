#include "atomic-density-analysis.h"

template <class T>
void atomic_density_analysis<T>::Setup (system_t& t) {

	AnalysisSet<system_t>::Setup(t);

	// grab the list of atomic names/types that will be used for the analysis and create the vector-histograms
	libconfig::Setting &atom_names = WaterSystem<T>::SystemParameterLookup("analysis.density.atom-names");
	for (int i = 0; i < atom_names.getLength(); i++)
	{
		std::string atom_name = atom_names[i];
		atom_name_list.push_back(atom_name);

		histogram_set hs (3, histogram_t(-10.0, WaterSystem<T>::posmax+10.0, 0.1));
		histograms.insert(histogram_map_elmt(atom_name, hs));
	}

	// narrow down the system atoms to just those with names we're looking for
	md_name_utilities::KeepByNames (t.int_atoms, atom_name_list);

}	// Setup


template <class T>
void atomic_density_analysis<T>::Analysis (system_t& t) { 

	//for (Atom_it it = t.int_atoms.begin(); it != t.int_atoms.end(); it++) {
		//VecR pos = (*it)->Position();
		//histograms[(*it)->Name()].push_back (pos);
	//}
		
	std::for_each (t.int_atoms.begin(), t.int_atoms.end(), std::bind2nd(binner, &histograms));
	
}

template <class T>
void atomic_density_analysis<T>::DataOutput (system_t& t) {

	rewind(t.Output());

	// first output the header of all the atom-names
	for (std::vector<std::string>::const_iterator it = atom_name_list.begin(); it != atom_name_list.end(); it++) {
		fprintf (t.Output(), "%s_x %s_y %s_z ", it->c_str(), it->c_str(), it->c_str());
	}
	fprintf (t.Output(), "\n");
	fflush(t.Output());

	// output the data from the histograms
	double dr = histograms.begin()->second[0].Resolution();
	double min = histograms.begin()->second[0].Min();
	double max = histograms.begin()->second[0].Max();

	for (double r = min; r < max; r+=dr) {
		for (std::vector<std::string>::const_iterator name = atom_name_list.begin(); name != atom_name_list.end(); name++) {

			histogram_set& hs = histograms[*name];
			for (int ax = 0; ax < 3; ax++) {
				fprintf (t.Output(), "% 8.3f % 8.3f ", r, hs[ax].Population(r)/t.timestep);
			}
		}
		fprintf (t.Output(), "\n");
	}

	fflush(t.Output());

}	// Data Output


// used to calculate the SFG spectrum based on the morita/hynes 2008 method
int main () {

	Analyzer<AmberSystem> analyzer;
	atomic_density_analysis<AmberSystem> analysis;
	analyzer.SystemAnalysis(analysis);

	return 0;
}










/*
	 template <class T, class U>
	 void DensityBinner<T,U>::operator() (T t) {
// if procesing a new atom type that hasn't yet been encountered, create a new histogram for it
if (!MapKeyExists(t))
AddNewHistogram(t);

_histograms[t->Name()][Analyzer<U>::PositionBin(t)]++;

return;
}

template <class T, class U>
// Output data from the histograms
void DensityBinner<T,U>::Output (FILE * output, const int timestep) {

rewind (output);

double scale = (!timestep) ? 1.0 : double(timestep);

// first print out a header consisting of column headers
fprintf (output, "Position");	// The position column
// column headers for each atom type (name)
for (Histogram_it it = _histograms.begin(); it != _histograms.end(); it++) {
fprintf (output, "%13s", (*it).first.c_str());
}
fprintf(output, "\n");

double position;
// for every position in the system
for (int pos = 0; pos < Analyzer<U>::posbins; pos++) {

// print the slab position
position = double(pos)*Analyzer<U>::posres + Analyzer<U>::posmin;
fprintf (output, "% 13.5f", position);

// go through each atom type
for (Histogram_it it = _histograms.begin(); it != _histograms.end(); it++) {
// and print the density value for the given position in it's own column
double density = double((*it).second[pos]) / scale;
fprintf (output, "% 13.5f", density);
}
fprintf(output, "\n");
}
return;
}

int main () {

libconfig::Config cfg;
cfg.readFile("system.cfg");

std::string filename = cfg.lookup("analysis.density.filename");
libconfig::Setting &analysis = cfg.lookup("analysis");
analysis.add("filename", libconfig::Setting::TypeString) = filename;

WaterSystemParams wsp (cfg);

DensityAnalyzer<GMXSystem> da (wsp);

da.SystemAnalysis();

return 0;
}
*/
