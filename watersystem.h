#ifndef	WATERSYSTEM_H_
#define	WATERSYSTEM_H_

#include "mdsystem.h"
#include "ambersystem.h"
#include "xyzsystem.h"
#include "utility.h"
#include "graph.h"

// Predicates (defined below)
template <class T> bool Name_pred (const T t, const std::string name);
template <class T> bool AtomicPosition_pred (const T t, const std::pair<double,double> extents);

struct WaterSystemParams {
  WaterSystemParams (		/* Constructor and initial/default values for system parameters */
		     std::string _output = "temp.dat",
		     const int _timesteps = 200000,
		     const bool _avg = false,
		     const coord _axis = y,
		     const VecR _ref_axis = VecR (0.0, 1.0, 0.0),
		     const int _output_freq = 100, const int _restart = 0,
		     const double _posmin = -20.0, const double _posmax = 150.0, 
		     const double _posres = 0.1,
		     const double _angmin = -1.0, const double _angmax = 1.0, 
		     const double _angres = 0.01,
		     const double _pbcflip = 20.0
		     ) :	/* The initialization of the system parameters */
  avg(_avg), output_filename(_output), output(fopen(_output.c_str(), "w")),
    axis(_axis), ref_axis(_ref_axis), output_freq(_output_freq), timesteps(_timesteps), restart(_restart),
    posmin(_posmin), posmax(_posmax), posres(_posres), 
    posbins (int((posmax-posmin)/posres)),
    pbcflip(_pbcflip),
    angmin(_angmin), angmax(_angmax), angres(_angres),
    angbins (int((angmax-angmin)/angres))
  { }

  bool avg;			/* Will averaging of two interfaces be performed? Can also be used for other functionality */

  std::string output_filename;
  FILE * output;		/* The file name and stream of the data output */

  coord axis;			/* The reference axis (generally normal to the interface */
  VecR ref_axis;

  int output_freq;		/* How often data is output to disk and info is posted to the screen */
  int timesteps, restart;	/* Number of timesteps to process, and restart = number of timesteps to skip */

  double posmin, posmax, posres; /* For generating histograms position data histograms */
  int posbins;			 /* The number of bins in the position histograms */
  double pbcflip;
  double angmin, angmax, angres; /* For generating histograms to bin angle data (min, max, bin width/resolution) */
  int angbins;

  void PrintParameters () {
	printf ("output_file = %s\naxis = %d, ref_axis = ", output_filename.c_str(), axis);
	ref_axis.Print();
	printf ("output_freq = %d, timesteps = %d, restart = %d\n", output_freq, timesteps, restart);
	printf ("posmin = % 8.3f, posmax = % 8.3f, posres = % 8.3f, posbins = %d\n", posmin, posmax, posres, posbins);
	printf ("pbcflip = % 8.3f\nangmin = % 8.3f, angmax = % 8.3f, angres = % 8.3f, angbins = %d\n",
		pbcflip, angmin, angmax, angres, angbins);
  }
};


template<class T>
class WaterSystem {

 public:

  WaterSystem (const WaterSystemParams& params);
  //WaterSystem (const int argc, const char **argv, const WaterSystemParams& params);
  virtual ~WaterSystem ();

  static WaterSystemParams wsp;

  static double	posmin, posmax;
  static double	pbcflip;			// location to flip about periodic boundaries
  static coord axis;					// axis normal to the infterface
  static double int_low, int_high, middle;		// the positions of analysis cutoffs

  static Water_ptr_vec	int_wats;		// interfacial waters, or just all the waters in the system depending on the function call
  static Mol_ptr_vec 	int_mols;
  static Atom_ptr_vec	int_atoms;		// interfacial water atoms (or as above)


  void OpenFile ();
  void Debug (string msg) const;

  void LoadAll ();							// Loads all the molecules and atoms in the system into the containers
  void FindWaters ();						// Find all the waters
  void FindMols (const string name);		// find all molecules with this name
  void FlipWaters (const coord axis = y);
  void SliceAtoms (const double low, const double high);
  void SliceWaters (const double low, const double high);
  void SliceWaterCoordination (const coordination coord);
  void FindInterfacialWaters ();

  void UpdateGraph () { _graph.UpdateGraph (int_atoms); }

 protected:

  T * sys;

  BondGraph	_graph;
};

template<typename T> WaterSystemParams WaterSystem<T>::wsp;

template<typename T> double WaterSystem<T>::posmin;
template<typename T> double WaterSystem<T>::posmax;
template<typename T> double WaterSystem<T>::pbcflip;
template<typename T> coord WaterSystem<T>::axis;

template<typename T> Water_ptr_vec WaterSystem<T>::int_wats;
template<typename T> Mol_ptr_vec WaterSystem<T>::int_mols;
template<typename T> Atom_ptr_vec WaterSystem<T>::int_atoms;

template <class T>
WaterSystem<T>::WaterSystem (const WaterSystemParams& params)
{

  WaterSystem::wsp = params;

  WaterSystem::posmin = params.posmin;
  WaterSystem::posmax = params.posmax;
  WaterSystem::pbcflip = params.pbcflip;
  WaterSystem::axis = params.axis;

  return;
}

template <class T>
WaterSystem<T>::~WaterSystem () {

  return;
}

template <class T>
void WaterSystem<T>::FindInterfacialWaters () {

  int_wats.clear();
  int_atoms.clear();

  Molecule * pmol;

  // go through the system
  for (int i = 0; i < sys->NumMols(); i++) {
    // grab each molecule
    pmol = sys->Molecules(i);

    // we're only looking at waters for SFG analysis right now
    if (pmol->Name() != "h2o") continue;

    // first thing we do is grab a water molecule to work with
    Water * water = static_cast<Water *>(pmol);

    double position = water->Atoms(0)->Position()[axis];
    // and find molecules that sit within the interface.
    if (position < pbcflip) position += MDSystem::Dimensions()[axis];	// adjust for funky boundaries
    // these values have to be adjusted for each system
    if (position < posmin or position > posmax) continue;				// set the interface cutoffs

    int_wats.push_back (water);
    RUN2(water->Atoms()) {
      int_atoms.push_back (water->Atoms(j));
    }
  }

  return;
}

template <class T>
void WaterSystem<T>::FlipWaters (const coord axis) {

  RUN (int_wats) {
    int_wats[i]->Flip(axis);
  }

  return;
}

template <class T>
void WaterSystem<T>::FindMols (const string name) {

  int_mols.clear();
  int_atoms.clear();

  Molecule * pmol;

  // go through the system
  for (int i = 0; i < sys->NumMols(); i++) {
    // grab each molecule
    pmol = sys->Molecules(i);

    // we're only looking at waters for SFG analysis right now
    if (pmol->Name() != name) continue;

    // add it to the group
    int_mols.push_back (pmol);

    // and then add all of its atoms
    RUN2(pmol->Atoms()) {
      int_atoms.push_back (pmol->Atoms(j));
    }
  }

  if ((int)int_mols.size() == 0) {
    printf ("\n\n***  Couldn't find any molecules with the residue name \"%s\" in the system\nCheck the system settings\n", name.c_str());
    exit(1);
  }

  fflush(stdout);
  return;
}

template <class T>
void WaterSystem<T>::LoadAll () {

  int_mols.clear();
  int_atoms.clear();

  Molecule * pmol;

  // go through the system
  for (int i = 0; i < sys->NumMols(); i++) {
    // grab each molecule to be added
	pmol = sys->Molecules(i);
    int_mols.push_back (pmol);
    // and then add all of its atoms
    RUN2(pmol->Atoms()) {
      int_atoms.push_back (pmol->Atoms(j));
    }
  }

  return;
}

template <class T>
void WaterSystem<T>::FindWaters () {

  int_wats.clear();
  int_atoms.clear();

  Molecule * pmol;
  Water * water;

  // go through the system
  for (int i = 0; i < sys->NumMols(); i++) {
    // grab each molecule
    pmol = sys->Molecules(i);

    // we're only looking at waters for SFG analysis right now
    if (pmol->Name() != "h2o") continue;

    // first thing we do is grab a water molecule to work with
    water = static_cast<Water *>(pmol);
    // add it to the group
    int_wats.push_back (water);

    // and then add all of its atoms
    RUN2(water->Atoms()) {
      int_atoms.push_back (water->Atoms(j));
    }
  }

  return;
}


// Let's the analysis only look at a particular piece of the system instead of the entire system. That is, only use waters that are between certain positions on the long-axis
template <class T>
void WaterSystem<T>::SliceWaters (const double low, const double high) {

  std::vector<Water *> wats;
  RUN (int_wats) {
    Water * wat = int_wats[i];
    Atom * oxy = wat->GetAtom ("O");
    VecR r = oxy->Position();
    double position = r[axis];
    if (position < pbcflip) {
      position += MDSystem::Dimensions()[axis];		// deal with the periodic cutoffs
    }

    if (position > low && position < high) {
      wats.push_back(wat);
    }
  }

  int_wats.clear();
  int_atoms.clear();

  RUN (wats) {
    int_wats.push_back(wats[i]);
    RUN2(wats[i]->Atoms()) {
      int_atoms.push_back(wats[i]->Atoms(j));
    }
  }

  return;
}

template <class T>
void WaterSystem<T>::SliceWaterCoordination (const coordination coord) {

  Water_ptr_vec wats;

  RUN (int_wats) {
    Water * wat = int_wats[i];
    coordination c = _graph.WaterCoordination(wat);
    if (c == coord) {
      wats.push_back(wat);
    }
  }

  int_wats.clear();
  int_atoms.clear();

  RUN (wats) {
    int_wats.push_back(wats[i]);
    RUN2(wats[i]->Atoms()) {
      int_atoms.push_back(wats[i]->Atoms(j));
    }
  }

  return;
}

/* used as a predicate for names equivalency */
template <class T>
bool Name_pred (const T t, const std::string name)
{
  return t->Name() == name;
}

/* predicate for finding a name in a list */
template <class T>
bool Name_in_list_pred (T t, std::vector<string> names) 
{
  return names.end() != find(names.begin(), names.end(), t->Name());
}

/* a predicate for name sorting */
template <class T>
struct Name_sort_pred {
    bool operator()(const T &left, const T &right) {
        return left->Name() < right->Name();
    }
};

/* removes all the atoms/molecules that are in the list of names to remove */
#define REMOVE_BY_NAMES(vec,type,names)	\
  std::sort(vec.begin(), vec.end(), Name_sort_pred<Atom *>());	\
  D_REMOVE_IF(vec, std::bind2nd(std::pointer_to_binary_function<type, std::vector<string>, bool>(Name_in_list_pred), names))

/* keeps all the things in vec with names that appear in the names-list. names should be a std::vector */
#define KEEP_BY_NAMES(vec,type,names)	\
  std::sort(vec.begin(), vec.end(), Name_sort_pred<Atom *>());	\
  D_REMOVE_IF(vec, std::bind2nd(not2(std::pointer_to_binary_function<type, std::vector<string>, bool>(Name_in_list_pred)), names))

// Removes any of the members with the given name
#define REMOVE_BY_NAME(vec,type,name)	\
  D_REMOVE_IF(vec,std::bind2nd(pointer_to_binary_function<const type, const std::string, bool>(Name_pred), name))

// a remove-if-not routine to remove anything *without* the given name
#define KEEP_BY_NAME(vec,type,name)	\
  D_REMOVE_IF(vec, std::bind2nd(not2(pointer_to_binary_function<const type, const std::string, bool>(Name_pred)), name))

/* used as a predicate to determine if a given atom is within the bounds specified */
template <class T>
bool AtomicPosition_pred (const T t, const std::pair<double,double> extents)
{
  double position = t->Position()[WaterSystem<AmberSystem>::axis];
  return position > extents.first && position < extents.second;
}

// Removes all elements not within the given slice bounded by the min and max positions within the watersystem slab
#define SLICE_BY_POSITION(vec,type,min,max)	\
  D_REMOVE_IF(vec, std::bind2nd(std::not2(std::pointer_to_binary_function<const type, const std::pair<double,double>, bool>(AtomicPosition_pred)), make_pair(min,max)))

#endif
