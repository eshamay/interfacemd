#ifndef TEST_H_
#define TEST_H_

#define EIGEN_MATRIXBASE_PLUGIN "/home/eshamay/md/src/EigenMatrixAddon.h"
#include <Eigen/Core>
USING_PART_OF_NAMESPACE_EIGEN

#include "../analysis.h"
#include "../utility.h"


class Tester : public Analyzer<XYZSystem>
{

  public:
	Tester (WaterSystemParams& params);

  private:

	void Setup ();
	void Analysis ();
	void DataOutput (const unsigned int timestep);
	void PostAnalysis () { return; }

	typedef std::pair<double, MolPtr>	mol_distance;
	std::vector< mol_distance > min_distances;
};

#endif

