#ifndef	WATERSYSTEM_H_
#define	WATERSYSTEM_H_

#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include "mdsystem.h"
#include "ambersystem.h"

#ifdef XYZ_SYS
#include "xyzsystem.h"
#endif

#ifdef GROMACS_SYS
#include "gmxsystem.h"
#endif


//#include "utility.h"
#include "graph.h"
#include <libconfig.h++>
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include <functional>


// easy way to carry around lots of config-file parameters
struct WaterSystemParams {
  WaterSystemParams () { }

  WaterSystemParams (libconfig::Config::Config * cfg) 
	: 
	  config_file (cfg), mpisys(false)
  { _Initialize(); }


  WaterSystemParams (libconfig::Config::Config * cfg, const bool mpisystem)
	: 
	  config_file (cfg), mpisys(true)
  { _Initialize(); }


  void _Initialize ()
  {
	try {
	  posmin = (config_file->lookup("analysis.position-range")[0]);
	  posmax = (config_file->lookup("analysis.position-range")[1]);

	  avg = config_file->lookup("analysis.averaging");

	  output_filename = (const char *)config_file->lookup("analysis.filename");
	  if (!mpisys || (mpisys && !mpi_world.rank())) {
		output = fopen(output_filename.c_str(), "w");
	  }

	  axis = (coord)((int)config_file->lookup("analysis.reference-axis"));
	  ref_axis = VecR(
		  config_file->lookup("analysis.reference-vector")[0], 
		  config_file->lookup("analysis.reference-vector")[1], 
		  config_file->lookup("analysis.reference-vector")[2]);
	  output_freq = config_file->lookup("analysis.output-frequency");
	  timesteps = config_file->lookup("system.timesteps");
	  restart = config_file->lookup("analysis.restart-time");
	  posres = config_file->lookup("analysis.resolution.position");
	  posbins  = int((posmax-posmin)/posres);
	  pbcflip = config_file->lookup("analysis.PBC-flip");
	  angmin = config_file->lookup("analysis.angle-range")[0];
	  angmax = config_file->lookup("analysis.angle-range")[1];
	  angres = config_file->lookup("analysis.resolution.angle");
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
  } // Initialize

  libconfig::Config * config_file;	/* Configuration file */
  boost::mpi::communicator	mpi_world;
  bool mpisys;

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

	WaterSystemParams wsp;
	static BondGraph graph;

	double	posmin, posmax;
	double	pbcflip;			// location to flip about periodic boundaries
	coord axis;				// axis normal to the infterface
	double int_low, int_high, middle;	// the positions of analysis cutoffs

	Atom_ptr_vec	sys_atoms;		// all atoms/mols in the system
	Mol_ptr_vec	sys_mols;
	Mol_ptr_vec	int_wats;		// interfacial waters, or just all the waters in the system depending on the function call
	Mol_ptr_vec 	int_mols;
	Atom_ptr_vec	int_atoms;		// interfacial water atoms (or as above)

	void OpenFile ();
	void Debug (string msg) const;

	bool MPI_Master () const { return !this->wsp.mpi_world.rank(); }
	int MPI_Rank () const { return this->wsp.mpi_world.rank(); }
	int MPI_Size () const { return this->wsp.mpi_world.size(); }
	void MPI_Print () const { cout << this->MPI_Rank() << "/" << this->MPI_Size() << endl; }

	virtual bool MPIStatus () const 
	{ return !this->wsp.mpisys || (this->wsp.mpisys && this->MPI_Master()); }




	static double AxisPosition (Atom * a) {
	  double pos = a->Position()[axis];
	  pos = (pos > pbcflip) ? pos : pos + MDSystem::Dimensions()[axis];
	  return pos;
	}

	typedef std::pair<double,double> Double_pair;
	// quick way to make a pair for the oft-used extents std::pair defaulting to the posmin/posmax in the config file
	/*
	   Double_pair ExtentPair (
	   const double low = this->posmin,
	   const double high = this->posmax) const {
	   return std::make_pair<double,double> (low, high);
	   }
	 */

	void LoadAll ();							// Loads all the molecules and atoms in the system into the containers
	void SliceWaterCoordination (const BondGraph::coordination c);

	// predicate determines if an atom sits within a particular slice of the system
	class AtomPositionInSlice : public std::binary_function<Atom *, Double_pair&, bool> {
	  public:
		bool operator() (Atom * atom, const Double_pair& extents) const
		{
		  double pos = AxisPosition (atom);
		  return pos > extents.first && pos < extents.second;
		}
	};

	// remove all the atoms in the defined slice (extents=(min,max)) of the slab
	void RemoveAtomsInSlice (Atom_ptr_vec& atoms, Double_pair& extents) {
	  atoms.erase(
		  remove_if(atoms.begin(), atoms.end(), std::bind2nd(AtomPositionInSlice(), extents)), atoms.end());
	  return;
	}

	void KeepAtomsInSlice (Atom_ptr_vec& atoms, Double_pair& extents) {
	  atoms.erase(
		  remove_if(atoms.begin(), atoms.end(), std::not1(std::bind2nd(AtomPositionInSlice(), extents))), atoms.end());
	  return;
	}

	// predicate to find if a water molecule is within a region (based on the position of the oxygen atom)
	class WaterInSlice : public std::binary_function<Molecule *, Double_pair, bool> {
	  private:
		AtomPositionInSlice apis;
	  public:
		bool operator() (const Molecule * wat, const Double_pair& extents) const
		{
		  return apis (wat->GetAtom("O"), extents);
		}
	};

	// This slices the water molecules and leaves only those within a given location of the slab. However, it does not do any loading or unloading of the atoms from int_atoms or anything else...
	void SliceWaters (Mol_ptr_vec& mols, Double_pair& extents) {
	  mols.erase(
		  remove_if(mols.begin(), mols.end(), std::not1(std::bind2nd(WaterInSlice(), extents))), mols.end());

	  this->UpdateAtoms(int_wats, int_atoms);
	  return;
	}

	// predicate to determine if a molecule is a water
	class IsWater_p : public std::unary_function<Molecule *, bool> {
	  public:
		bool operator() (const Molecule * mol) const
		{ return mol->Name() == "h2o"; }
	};

	/* loads the int_wats and int_atoms with only waters and water atoms */
	void LoadWaters () {
	  LoadAll();
	  int_wats.clear();

	  // copy over all the water molecules into the int_wats container
	  std::remove_copy_if(		// have to use remove_copy_if because the STL doesn't have a copy_if!!!
		  sys_mols.begin(), sys_mols.end(), std::back_inserter(int_wats), std::not1(IsWater_p()));

	  // load in the water atoms to int_atoms
	  this->UpdateAtoms (int_wats, int_atoms);

	  return;
	}

	// predicate to test if the name of an atom or molecule (determined by the template parameter) is found in a vector of names
	template <class U>
	  class NameInList : public std::binary_function<U, std::vector<string>,bool> {
		public:
		  bool operator() (const U u, const std::vector<string>& names) const
		  {
			return names.end() != std::find(names.begin(), names.end(), u->Name());
		  }
	  };

	// predicate to determine if the given atom or molecule has the given name
	template <class U>
	  class IsName : public std::binary_function<U, std::string, bool> {
		public:
		  bool operator() (const U u, const std::string name) const
		  {
			return u->Name() == name;
		  }
	  };

	template <class U>
	  void KeepByName (U& u, std::string& name) {
		u.erase(
			remove_if(u.begin(), u.end(), not1(std::bind2nd(IsName<typename U::value_type>(), name))), u.end()
			);
		return;
	  }

	// keep elements of the vector with names matching one of those in the list
	template <class U>
	  void KeepByNames (U& u, std::vector<string>& names) {
		u.erase(
			remove_if(u.begin(), u.end(), not1(std::bind2nd(NameInList<typename U::value_type>(), names))), u.end());
		return;
	  }


	// remove all elements that have the given name
	template <class U>
	  void RemoveByName (U& u, std::string& name) {
		u.erase(
			remove_if(u.begin(), u.end(), std::bind2nd(IsName<typename U::value_type>(), name)), u.end()
			);
		return;
	  }

	// remove all elements that have names matching any of those in the list of names supplied
	template <class U>
	  void RemoveByNames (U& u, std::vector<string>& names) {
		u.erase(
			remove_if(u.begin(), u.end(), std::bind2nd(NameInList<typename U::value_type>(), names)), u.end());
		return;
	  }

	// Predicate to test if a water molecule has a given coordination (H-bonding pattern)
	class WaterCoordination_p : public std::binary_function<Water *, BondGraph::coordination, bool> {
	  public:
		bool operator() (const Water * wat, const BondGraph::coordination c) const {
		  return graph.WaterCoordination(wat) == c;
		}
	};

	void KeepWaterByCoordination (Mol_ptr_vec& mols, const BondGraph::coordination c) {
	  mols.erase(
		  remove_if(mols.begin(), mols.end(), std::not1(std::bind2nd(WaterCoordination_p(), c))), mols.end());
	  return;
	}

	// Sets the atom container to hold only the atoms of the set of molecules given
	void UpdateAtoms (const Mol_ptr_vec& mols, Atom_ptr_vec& atoms) {
	  atoms.clear();
	  for (Mol_it it = mols.begin(); it != mols.end(); it++) {
		std::copy ((*it)->begin(), (*it)->end(), std::back_inserter(atoms));
	  }
	  return;
	}

	void UpdateGraph () { graph.UpdateGraph (int_atoms); }

  protected:

	T * sys;

};

template<typename T> BondGraph WaterSystem<T>::graph;

  template <class T>
WaterSystem<T>::WaterSystem (const WaterSystemParams& params)
  :
	posmin(params.posmin), posmax(params.posmax), pbcflip(params.pbcflip), axis(params.axis), wsp(params)
{ return; }


template <class T>
WaterSystem<T>::~WaterSystem () { return; }



template <class T>
void WaterSystem<T>::LoadAll () {

  sys_mols.clear();
  sys_atoms.clear();
  int_mols.clear();
  int_atoms.clear();

  Mol_ptr_vec mols = sys->Molecules();
  // copy the molecules and atoms into each container
  std::copy(mols.begin(), mols.end(), std::back_inserter(sys_mols));
  std::copy(mols.begin(), mols.end(), std::back_inserter(int_mols));

  // copy the atoms into the atom containers
  this->UpdateAtoms (sys_mols, sys_atoms);
  this->UpdateAtoms (int_mols, int_atoms);

  return;
}

#endif
