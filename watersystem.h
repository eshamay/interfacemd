#ifndef	WATERSYSTEM_H_
#define	WATERSYSTEM_H_

#include "mdsystem.h"
#include "ambersystem.h"
#include "xyzsystem.h"
#include "utility.h"
#include "graph.h"
#include <libconfig.h++>
#include <cstdlib>
#include <iomanip>

// Predicates (defined below)
template <class T> bool Name_pred (const T t, const std::string name);
template <class T> bool AtomicPosition_pred (const T t, const std::pair<double,double> extents);

struct WaterSystemParams {
  /*
  WaterSystemParams (		// Constructor and initial/default values for system parameters
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
		     ) :	// The initialization of the system parameters
  avg(_avg), output_filename(_output), output(fopen(_output.c_str(), "w")),
    axis(_axis), ref_axis(_ref_axis), output_freq(_output_freq), timesteps(_timesteps), restart(_restart),
    posmin(_posmin), posmax(_posmax), posres(_posres), 
    posbins (int((posmax-posmin)/posres)),
    pbcflip(_pbcflip),
    angmin(_angmin), angmax(_angmax), angres(_angres),
    angbins (int((angmax-angmin)/angres))
  { }
  */

  WaterSystemParams () { }

  WaterSystemParams (libconfig::Config::Config& cfg)
  {
    try {
      config_file = &cfg;
      posmin = (cfg.lookup("analysis.position-range")[0]);
      posmax = (cfg.lookup("analysis.position-range")[1]);

      avg = cfg.lookup("analysis.averaging");
      output_filename = (const char *)cfg.lookup("analysis.filename");
      output = fopen(output_filename.c_str(), "w");
      axis = (coord)((int)cfg.lookup("analysis.reference-axis"));
      ref_axis = VecR(
	  cfg.lookup("analysis.reference-vector")[0], 
	  cfg.lookup("analysis.reference-vector")[1], 
	  cfg.lookup("analysis.reference-vector")[2]);
      output_freq = cfg.lookup("analysis.output-frequency");
      timesteps = cfg.lookup("system.timesteps");
      restart = cfg.lookup("analysis.restart-time");
      posres = cfg.lookup("analysis.resolution.position");
      posbins  = int((posmax-posmin)/posres);
      pbcflip = cfg.lookup("analysis.PBC-flip");
      angmin = cfg.lookup("analysis.angle-range")[0];
      angmax = cfg.lookup("analysis.angle-range")[1];
      angres = cfg.lookup("analysis.resolution.angle");
      angbins  = int((angmax-angmin)/angres);
    }
    catch(const libconfig::SettingTypeException &stex) {
      std::cerr << "Something is wrong with the configuration parameters or file - check syntax\n(watersystem.h)" << std::endl;
      exit(EXIT_FAILURE);
    }
    catch(const libconfig::SettingNotFoundException &snfex) {
      std::cerr << "A setting is missing from the configuration file!" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  libconfig::Config * config_file;	/* Configuration file */

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
  static coord axis;				// axis normal to the infterface
  static double int_low, int_high, middle;	// the positions of analysis cutoffs

  static Atom_ptr_vec	sys_atoms;		// all atoms/mols in the system
  static Mol_ptr_vec	sys_mols;
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
  void SliceWaterCoordination (const BondGraph::coordination c);
  void FindInterfacialWaters ();

  void KeepAtomsByNames (Atom_ptr_vec& atoms, std::vector<string>& names);
  void RemoveAtomsByNames (Atom_ptr_vec& atoms, std::vector<string>& names);

  void RemoveAtomsInSlice (Atom_ptr_vec& atoms, std::pair<double,double>& extents);
  void KeepAtomsInSlice (Atom_ptr_vec& atoms, std::pair<double,double>& extents);

  void UpdateGraph () { _graph.UpdateGraph (int_atoms); }

  void PrintAtoms (const Atom_ptr_vec& atoms) const
  {
    for (Atom_it it = atoms.begin(); it != atoms.end(); it++) 
    { (*it)->Print(); }
  }

 protected:

  T * sys;

  BondGraph	_graph;

  class AtomPositionInSlice : public std::binary_function<Atom *,std::pair<double,double>,bool> {
    public:
      bool operator() (const Atom * atom, const std::pair<double,double>& extents) const
      {
	double pos = atom->Position()[axis];
	return pos > extents.first && pos < extents.second;
      }
  };

  class AtomNameInList : public std::binary_function<Atom *,std::vector<string>,bool> {
    public:
      bool operator() (const Atom * atom, const std::vector<string>& names) const
      {
	return names.end() != std::find(names.begin(), names.end(), atom->Name());
      }
  };




};

template<typename T> WaterSystemParams WaterSystem<T>::wsp;

template<typename T> double WaterSystem<T>::posmin;
template<typename T> double WaterSystem<T>::posmax;
template<typename T> double WaterSystem<T>::pbcflip;
template<typename T> coord WaterSystem<T>::axis;

template<typename T> Atom_ptr_vec WaterSystem<T>::sys_atoms;
template<typename T> Mol_ptr_vec WaterSystem<T>::sys_mols;
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

  sys_mols.clear();
  sys_atoms.clear();
  int_mols.clear();
  int_atoms.clear();

  Molecule * pmol;

  // go through the system
  for (int i = 0; i < sys->NumMols(); i++) {
    // grab each molecule to be added
    pmol = sys->Molecules(i);
    int_mols.push_back (pmol);
    sys_mols.push_back (pmol);
    // and then add all of its atoms
    RUN2(pmol->Atoms()) {
      int_atoms.push_back (pmol->Atoms(j));
      sys_atoms.push_back (pmol->Atoms(j));
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
void WaterSystem<T>::SliceWaterCoordination (const BondGraph::coordination c) {

  Water_ptr_vec wats;

  RUN (int_wats) {
    Water * wat = int_wats[i];
    BondGraph::coordination cd = _graph.WaterCoordination(wat);
    if (cd == c) {
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

/* keep only atoms in the given atom_ptr_vec that have a name in the names list */
template <class T>
void WaterSystem<T>::KeepAtomsByNames (Atom_ptr_vec& atoms, std::vector<string>& names) {
  atoms.erase(
      remove_if(atoms.begin(), atoms.end(), 
	not1(std::bind2nd(AtomNameInList(), names))), 
      atoms.end());
  return;
}

/* remove all atoms that have a name in the list */
template <class T>
void WaterSystem<T>::RemoveAtomsByNames (Atom_ptr_vec& atoms, std::vector<string>& names) {
  atoms.erase(
      remove_if(atoms.begin(), atoms.end(),
	std::bind2nd(AtomNameInList(), names)),
      atoms.end());
  return;
}

/* remove all the atoms in the defined slice (extents=(min,max)) of the slab */
template <class T>
void WaterSystem<T>::RemoveAtomsInSlice (Atom_ptr_vec& atoms, std::pair<double,double>& extents) {
  atoms.erase(
      remove_if(atoms.begin(), atoms.end(),
	std::bind2nd(AtomPositionInSlice(), extents)),
      atoms.end());
  return;
}

template <class T>
void WaterSystem<T>::KeepAtomsInSlice (Atom_ptr_vec& atoms, std::pair<double,double>& extents) {
  atoms.erase(
      remove_if(atoms.begin(), atoms.end(),
	std::not1(std::bind2nd(AtomPositionInSlice(), extents))),
      atoms.end());
  return;
}

#endif
