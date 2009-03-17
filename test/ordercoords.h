#ifndef	ORDER_PARAMS_H_
#define	ORDER_PARAMS_H_

#include "../watersystem.h"

#define AVG

typedef std::vector<double> 		Double_vec;
typedef std::vector<Double_vec> 	Double_histo;
typedef std::vector<int>			Int_vec;
typedef std::vector<Int_vec> 		Int_histo;

class CoordOrderParams : public WaterSystem {

public:

	CoordOrderParams (const int argc, const char **argv, const WaterSystemParams& params);

	void Analysis ();

protected:

	// working with order parameters
	Double_histo S1;
	Double_histo S2_num;
	Double_histo S2_den;
	Int_histo number_density;

	double	angmax;
	double	angmin;
	double	angres;
	int		angbins;

	// working with coordination types
	std::vector<coordination> vcoords;
	coord_map name_map;

	void SetHistograms (const int argc, const char **argv);
	void BinOrderParameter (Water * wat);
	void OutputData ();

};

#endif
