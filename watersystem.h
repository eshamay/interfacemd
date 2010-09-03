#pragma once
#ifndef	WATERSYSTEM_H_
#define	WATERSYSTEM_H_

#include "mdsystem.h"
#include "ambersystem.h"
#include "xyzsystem.h"

#ifdef GROMACS_SYS
#include "gmxsystem.h"
#endif

#include "utility.h"

#include <libconfig.h++>
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <numeric>


// easy way to carry around lots of config-file parameters
class WaterSystemParams {
	public: 
		WaterSystemParams () { }

		~WaterSystemParams () { 
			//std::cout << "WaterSystemParams dtor" << std::endl;
			delete config_file; 
		}

		WaterSystemParams (std::string configuration_filename)
		{
			try {

				config_file = new libconfig::Config();
				printf ("\nUsing configuration file: \"%s\"\n", configuration_filename.c_str());
				config_file->readFile(configuration_filename.c_str());

				posmin = (config_file->lookup("analysis.position-range")[0]);
				posmax = (config_file->lookup("analysis.position-range")[1]);

				avg = config_file->lookup("analysis.averaging");
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

		}

		libconfig::Config * config_file;	/* Configuration file */

		bool avg;			/* Will averaging of two interfaces be performed? Can also be used for other functionality */

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
			//printf ("output_file = %s\naxis = %d, ref_axis = ", output_filename.c_str(), axis);
			printf ("axis = %d, ref_axis = ", axis);
			ref_axis.Print();
			printf ("output_freq = %d, timesteps = %d, restart = %d\n", output_freq, timesteps, restart);
			printf ("posmin = % 8.3f, posmax = % 8.3f, posres = % 8.3f, posbins = %d\n", posmin, posmax, posres, posbins);
			printf ("pbcflip = % 8.3f\nangmin = % 8.3f, angmax = % 8.3f, angres = % 8.3f, angbins = %d\n",
					pbcflip, angmin, angmax, angres, angbins);
		}
};


class AmberAnalysisSet;
class XYZAnalysisSet;

template<class T>
class WaterSystem {

	public:

		WaterSystem (const std::string configuration_filename);
		virtual ~WaterSystem ();

		static WaterSystemParams * SystemParameters () { return wsp; }

		static libconfig::Setting& SystemParameterLookup (std::string param) {
			try {
				return WaterSystem<T>::wsp->config_file->lookup(param);
			}
			catch(const libconfig::SettingTypeException &stex) {
				std::cerr << "Something is wrong with the " << param << " setting in the system.cfg configuration file." << std::endl;
				exit(EXIT_FAILURE);
			}
			catch(const libconfig::SettingNotFoundException &snfex) {
				std::cerr << "Couldn't find the " << param << " setting in the system.cfg configuration file." << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		T * System () { return sys; }

		static double	posmin, posmax;
		static double	pbcflip;			// location to flip about periodic boundaries
		static coord axis;				// axis normal to the infterface
		static double int_low, int_high, middle;	// the positions of analysis cutoffs

		static Atom_ptr_vec	sys_atoms;		// all atoms/mols in the system
		static Mol_ptr_vec	sys_mols;

		static Atom_ptr_vec	int_atoms;		// interfacial water atoms (or as above)
		static Mol_ptr_vec 	int_mols;

		static Mol_ptr_vec	int_wats;		// interfacial waters, or just all the waters in the system depending on the function call

		void OpenFile ();


		static double AxisPosition (const AtomPtr a) {
			double pos = a->Position()[axis];
			pos = (pos > pbcflip) ? pos : pos + MDSystem::Dimensions()[axis];
			return pos;
		}

		typedef std::pair<double,double> Double_pair;
		// quick way to make a pair for the oft-used extents std::pair defaulting to the posmin/posmax in the config file
		Double_pair ExtentPair (
				const double low = WaterSystem<T>::posmin,
				const double high = WaterSystem<T>::posmax) const {
			return std::make_pair<double,double> (low, high);
		}

		void LoadAll ();							// Loads all the molecules and atoms in the system into the containers
		void SliceWaterCoordination (const bondgraph::coordination c);

		// predicate determines if an atom sits within a particular slice of the system
		class AtomPositionInSlice : public std::binary_function<AtomPtr, Double_pair&, bool> {
			public:
				bool operator() (const AtomPtr atom, const Double_pair& extents) const
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
		template <typename U>
			class WaterInSlice : public std::binary_function<U, Double_pair, bool> {
				private:
					AtomPositionInSlice apis;
				public:
					bool operator() (const U wat, const Double_pair& extents) const
					{
						return apis (wat->GetAtom("O"), extents);
					}
			};

		// This slices the water molecules and leaves only those within a given location of the slab. However, it does not do any loading or unloading of the atoms from int_atoms or anything else...
		template <typename U>
			static void SliceWaters (std::vector<U>& mols, Double_pair& extents) {
				mols.erase(
						remove_if(mols.begin(), mols.end(), std::not1(std::bind2nd(WaterInSlice<U>(), extents))), mols.end());

				//this->UpdateAtoms(int_wats, int_atoms);
				return;
			}


		/* loads the int_wats and int_atoms with only waters and water atoms */
		void LoadWaters () {
			LoadAll();
			int_wats.clear();

			// copy over all the water molecules into the int_wats container
			algorithm_extra::copy_if (sys_mols.begin(), sys_mols.end(), std::back_inserter(int_wats), member_functional::mem_fun_eq(&Molecule::MolType, Molecule::H2O));

			// load in the water atoms to int_atoms
			this->UpdateAtoms (int_wats, int_atoms);

			return;
		}


		// Predicate to test if a water molecule has a given coordination (H-bonding pattern)
		class WaterCoordination_p : public std::binary_function<Water *, bondgraph::coordination, bool> {
			public:
				bool operator() (const Water * wat, const bondgraph::coordination c) const {
					return sys->graph.WaterCoordination(wat) == c;
				}
		};


		void KeepWaterByCoordination (Mol_ptr_vec& mols, const bondgraph::coordination c) {
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

		void UpdateGraph () { sys->graph.UpdateGraph (int_atoms); }
		bondgraph::BondGraph& Graph () const { return sys->graph; }


	protected:

		static WaterSystemParams * wsp;

		T * sys;

		virtual void _InitializeSystem ();

};

template<typename T> WaterSystemParams * WaterSystem<T>::wsp;

template<typename T> double WaterSystem<T>::posmin;
template<typename T> double WaterSystem<T>::posmax;
template<typename T> double WaterSystem<T>::pbcflip;
template<typename T> coord WaterSystem<T>::axis;

template<typename T> Atom_ptr_vec WaterSystem<T>::sys_atoms;
template<typename T> Mol_ptr_vec WaterSystem<T>::sys_mols;
template<typename T> Mol_ptr_vec WaterSystem<T>::int_wats;
template<typename T> Mol_ptr_vec WaterSystem<T>::int_mols;
template<typename T> Atom_ptr_vec WaterSystem<T>::int_atoms;


	template <class T>
WaterSystem<T>::WaterSystem (const std::string configuration_filename) 
{

	wsp = new WaterSystemParams(configuration_filename);

	posmin = wsp->posmin;
	posmax = wsp->posmax;
	pbcflip = wsp->pbcflip;
	axis = wsp->axis;

	try {
		this->_InitializeSystem();
	}
	catch (const libconfig::SettingTypeException &stex) {
		std::cerr << "WaterSystem<T>::_InitializeSystem() -- Something wrong with initializing the system. Try checking filenames in the system.cfg" << std::endl;
		exit(EXIT_FAILURE);
	}

	return;
}

template <class T>
WaterSystem<T>::~WaterSystem () {
	delete wsp;
	return;
}

template <>
void WaterSystem<AmberSystem>::_InitializeSystem () {
	try {
		this->sys = new AmberSystem(
				wsp->config_file->lookup("system.files.prmtop"),
				wsp->config_file->lookup("system.files.mdcrd"),
				wsp->config_file->lookup("system.files.mdvel"));
	}
	catch (const libconfig::SettingNotFoundException &snfex) {
		std::cerr << "Couldn't find the Amber system filenames listed in the configuration file" << std::endl;
		exit(EXIT_FAILURE);
	}
	return;
}

template <>
void WaterSystem<XYZSystem>::_InitializeSystem () {

	try {
		std::string filepath = wsp->config_file->lookup("system.files.xyzfile");

		double a,b,c;
		a = wsp->config_file->lookup("system.dimensions")[0];
		b = wsp->config_file->lookup("system.dimensions")[1];
		c = wsp->config_file->lookup("system.dimensions")[2];

		std::string wanniers = wsp->config_file->lookup("system.files.wanniers");
		VecR dims(a,b,c);
		printf ("system dimensions are: ");
		dims.Print();

		this->sys = new XYZSystem(filepath, dims, wanniers);
	}
	catch (const libconfig::SettingNotFoundException &snfex) {
		std::cerr << "Couldn't find the xyz system parameters in the configuration file" << std::endl;
		exit(EXIT_FAILURE);
	}

	return;
}

#ifdef GROMACS_SYS
template <>
void WaterSystem<GMXSystem>::_InitializeSystem () {
	std::string trr = wsp->config_file->lookup("system.files.gmx-trrfile");
	std::string gro = wsp->config_file->lookup("system.files.gmx-grofile");
	this->sys = new GMXSystem(trr.c_str(), gro.c_str());
	return;
}
#endif


template <class T>
void WaterSystem<T>::LoadAll () {

	sys_mols.clear();
	sys_atoms.clear();
	int_mols.clear();
	int_atoms.clear();

	// copy the molecules and atoms into each container
	std::copy(sys->begin_mols(), sys->end_mols(), std::back_inserter(sys_mols));
	std::copy(sys->begin_mols(), sys->end_mols(), std::back_inserter(int_mols));

	std::copy(sys->begin(), sys->end(), std::back_inserter(sys_atoms));
	std::copy(sys->begin(), sys->end(), std::back_inserter(int_atoms));

	return;
}




#endif
