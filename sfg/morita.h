#ifndef MORITA_H_
#define MORITA_H_

#include "analysis.h"
#include "sfgunits.h"
#include "moritah2o.h"
#include "watersfg.h"

#include <Eigen/LU>
#include <iostream>
#include <mkl_blas.h>
#include <mkl_lapack.h>
#include <functional>

namespace morita {

	USING_PART_OF_NAMESPACE_EIGEN

	typedef std::vector< std::complex<double> > Complex_vec;

	class MoritaAnalysis : public AnalysisSet< Analyzer<AmberSystem> > {
		public:
			MoritaAnalysis () : 
				AnalysisSet< Analyzer<AmberSystem> > (
						std::string("Analysis of the waters in an MD system using the technique of Morita/Hynes (2000) with a system produced by the Amber package (or any compliant dataset (prmtop, mdcrd, mdvel)"),
						std::string("amber-morita2000.dat")) { }

			~MoritaAnalysis ();

			typedef Analyzer<AmberSystem> system_t;

			virtual void Analysis (system_t&);
			virtual void DataOutput (system_t&);

		protected:

			static SFGCalculator	sfg;
			static Complex_vec TimestepChi;		// chi spectrum of several different molecules collected over an entire timestep
			static Complex_vec TotalChi;		// Running total of all the data for several timesteps

			static unsigned long numMolsProcessed;
			static bool firstMol;
			static bool firstTimeStep;

			std::vector<MoritaH2O *> all_wats;
			std::vector<MoritaH2O *> analysis_wats;

			void SetupSystemWaters (system_t& t);
			void SelectAnalysisWaters ();


			template <typename Iter1, typename Iter2>
				// Build up running totals of the chi spectra
				static void CollectChi (Iter1 pBegin, Iter1 pEnd, Iter2 pBegin2)
				{
					typedef typename std::iterator_traits<Iter1>::value_type value_type;

					std::transform(pBegin, pEnd, pBegin2, pBegin, std::plus<value_type>());
				}

			//static void CollectChi (Complex_vec& newchi, Complex_vec& totalchi);

			void FlipWaters (const coord axis = y);

			class SFGProcessor : public std::unary_function<MoritaH2O *,bool> {
				private:
					Complex_vec Molecular_Beta;			// chi for a given molecule

				public:

					bool operator() (MoritaH2O * water) {

						Molecular_Beta.clear();

						// and then calculate the chi spectrum for the molecule 
						// 0,2,1 = SSP. S = X and Z axes, P = Y axis
						Molecular_Beta = MoritaAnalysis::sfg.Beta (*water, 0,2,1);

						MoritaAnalysis::numMolsProcessed++;

						// when starting a new timestep...
						if (MoritaAnalysis::firstMol) {
							MoritaAnalysis::TimestepChi.resize (Molecular_Beta.size(), std::complex<double> (0.0, 0.0));
							MoritaAnalysis::firstMol = false;
						}

						// for the very first timestep...
						if (MoritaAnalysis::firstTimeStep) {
							MoritaAnalysis::TotalChi.clear();
							MoritaAnalysis::TotalChi.resize (Molecular_Beta.size(), std::complex<double> (0.0, 0.0));
							MoritaAnalysis::firstTimeStep = false;
						}

						// perform the summation for averaging spectra over all the water molecules
						MoritaAnalysis::CollectChi (MoritaAnalysis::TimestepChi.begin(), MoritaAnalysis::TimestepChi.end(), Molecular_Beta.begin());

						return true;
					}
			};

			SFGProcessor SFGProcess;
	};


}	// namespace morita

#endif
