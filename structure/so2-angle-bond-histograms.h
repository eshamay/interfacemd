#ifndef TEST_H_
#define TEST_H_

#define EIGEN_MATRIXBASE_PLUGIN "/home/eshamay/md/src/EigenMatrixAddon.h"
#include <Eigen/Core>
USING_PART_OF_NAMESPACE_EIGEN

#include "analysis.h"
#include "utility.h"


class Tester : public Analyzer<XYZSystem>
{
  public:
	Tester (WaterSystemParams& params);

  private:
	void Setup ();
	void Analysis ();
	void DataOutput ();
	void PostAnalysis ();

	SulfurDioxide * so2;
	typedef std::vector<double> double_vec;
	typedef double_vec::const_iterator double_it;
	double_vec so, theta;
};

#endif

