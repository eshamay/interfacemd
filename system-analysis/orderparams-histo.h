#ifndef	ORDER_PARAMS_H_
#define	ORDER_PARAMS_H_

//#define AVG
//#define RESTART

#include "../watersystem.h"

class OrderParameters : public WaterSystem<AmberSystem> {

public:

	OrderParameters (int argc, const char **argv, const WaterSystemParams& params);

	void Analysis ();

protected:

	double	angmax;
	double	angmin;
	double	angres;
	int		angbins;

	vector< vector<int>	> _data;

	void OutputData ();
	//void OutputStatus () const;

};

#endif
