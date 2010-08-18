#ifndef TEST_H_
#define TEST_H_

#define EIGEN_MATRIXBASE_PLUGIN "/home/eshamay/md/src/EigenMatrixAddon.h"
#include <Eigen/Core>
USING_PART_OF_NAMESPACE_EIGEN

#include "../analysis.h"
#include "../utility.h"


typedef std::vector<double> double_vec;
typedef double_vec::const_iterator double_it;

class StructureAnalyzer : public Analyzer<XYZSystem>
{
  public:
	StructureAnalyzer () 
	  : Analyzer<XYZSystem> () 
	{
	  PromptForAnalysisFunction(); 
	}

  protected:

	// output a bit of text to stdout and ask the user for a choice as to which type of analysis to perform - then do it.
	void PromptForAnalysisFunction ();
	// The set of possible analyses to perform on a given system
	std::vector<XYZAnalysisSet *> analyses;
};







class angle_bond_histogram_analyzer : public XYZAnalysisSet {
  public:
	angle_bond_histogram_analyzer (std::string desc, std::string fn) : XYZAnalysisSet (desc,fn) { }
	virtual ~angle_bond_histogram_analyzer() { }
	void DataOutput (system_t& t);
  protected:
	double_vec bondlengths, angles;
};

class h2o_angle_bond_histogram_analyzer : public angle_bond_histogram_analyzer {
  public:
	virtual ~h2o_angle_bond_histogram_analyzer() { }
	h2o_angle_bond_histogram_analyzer () :
	  angle_bond_histogram_analyzer(
		  std::string ("H2O H-O-H angle and O-H bondlength histograms"),
		  std::string ("h2o-bond-angle-histograms.dat")) { }

	void Setup (system_t& t) {
	  XYZAnalysisSet::Setup(t);
	}

	void Analysis (system_t& t);
  protected:
	Water * wat;

};

class so2_angle_bond_histogram_analyzer : public angle_bond_histogram_analyzer {
  public:
	virtual ~so2_angle_bond_histogram_analyzer() { }
	so2_angle_bond_histogram_analyzer() :
	  angle_bond_histogram_analyzer (
		  std::string("SO2 O-S-O angle and S-O bondlength histograms"),
		  std::string("so2-bond-angle-histograms.dat")) { }

	void Setup (system_t& t) {
	  XYZAnalysisSet::Setup(t);
	}

	void Analysis (system_t& t);
  protected:
	SulfurDioxide * so2;
};





class so2_angle_bond_analyzer : public XYZAnalysisSet {
  public:
	virtual ~so2_angle_bond_analyzer () { }
	so2_angle_bond_analyzer () :
	  XYZAnalysisSet (
		  std::string("SO2 molecular angle and S-O bondlengths"),
		  std::string ("so2-angle+bonds.dat")) { }

	void Setup (system_t& t) {
	  XYZAnalysisSet::Setup(t);
	  fprintf (t.Output(), "so1 so2 theta\n");
	}

	void Analysis (system_t& t);

  protected:
	SulfurDioxide * so2;
	double angle;
};




class so2_closest_atoms_analyzer : public XYZAnalysisSet {
  public:
	so2_closest_atoms_analyzer (std::string desc, std::string fn) : XYZAnalysisSet (desc,fn) { }
	virtual ~so2_closest_atoms_analyzer () { }
  protected:
	SulfurDioxide * so2;
	bondgraph::distance_vec closest;
};

class so2_closest_H_analyzer : public so2_closest_atoms_analyzer {
  public:
	virtual ~so2_closest_H_analyzer () { }
	so2_closest_H_analyzer () :
	  so2_closest_atoms_analyzer(
		  std::string("SO2 closest hydrogen analysis (reports the distance to the 3 hydrogens clsoest to each of the SO2 oxygens)"),
		  std::string("so2-closest-Hs.dat")) { }

	void Setup (system_t& t) {
	  XYZAnalysisSet::Setup(t);
	  fprintf (t.Output(), "o11 o12 o13 o21 o22 o23\n");
	}
	void Analysis (system_t& t);
  protected:
	AtomPtr o1,o2;
};

class so2_closest_O_analyzer : public so2_closest_atoms_analyzer {
  public:
	virtual ~so2_closest_O_analyzer () { }
	so2_closest_O_analyzer () :
	  so2_closest_atoms_analyzer(
		  std::string ("SO2 closest oxygen analysis (reports the distance to the 4 oxygens closest to the SO2 sulfur)"),
		  std::string ("so2-closest-Os.dat")) { }

	void Setup (system_t& t) {
	  XYZAnalysisSet::Setup(t);
	}
	void Analysis (system_t& t);
  protected:
	AtomPtr s;
};




class so2_hbond_factor_analyzer : public XYZAnalysisSet {
  public:
	virtual ~so2_hbond_factor_analyzer () { }
	so2_hbond_factor_analyzer () :
	  XYZAnalysisSet (
		  std::string ("H-sharing factor - unitless factor for studying the hydrogen-bond character of neighboring H's"),
		  std::string ("so2-hbond-factors.dat")) { }

	void Setup (system_t& t) {
	  XYZAnalysisSet::Setup (t);
	  fprintf(t.Output(), "timestep q1 q2\n");
	}

	void Analysis (system_t& t);
  protected:
	SulfurDioxide * so2;
	AtomPtr o1,o2;
};


#endif

