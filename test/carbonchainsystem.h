#ifndef	CARBONCHAINSYSTEM_H_
#define	CARBONCHAINSYSTEM_H_

#include "../watersystem.h"
#include "../carbonchain.h"
#include "../decane.h"
#include "../utility.h"

class CarbonChainSystem : public WaterSystem<AmberSystem> {

	public:

		CarbonChainSystem 
			(
			 int argc, 
			 const char **argv, 
			 const WaterSystemParams& params
			);

		vector<int> OrientationHistogram 
			(
			 const WaterSystemParams& params,
			 const vector<Molecule *> mols,
			 const string name,
			 VecR (*axisFunc)(const Atom * atom)
			)

		void Orientation-Analysis ();

};

#endif
