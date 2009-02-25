#ifndef	ORDER_PARAMS_H_
#define	ORDER_PARAMS_H_

#include "../watersystem.h"

#define AVG

class OrderParameters : public WaterSystem {

public:

	OrderParameters (int argc, char **argv, const WaterSystemParams& params);

	void Analysis ();

protected:

	vector<double> S1;
	vector<double> S2_num;
	vector<double> S2_den;
	vector<int> number_density;

	double	angmax;
	double	angmin;
	double	angres;
	int		angbins;

	void OutputData ();
	//void OutputStatus () const;

};

#endif
