#ifndef	ORDER_PARAMS_H_
#define	ORDER_PARAMS_H_

#define AVG
#define RESTART

#include "../watersystem.h"


class OrderParameters : public WaterSystem {

public:

	OrderParameters (int argc, const char **argv, const WaterSystemParams& params);

	void Analysis ();

protected:

/*
	vector<double> phi;
	vector<double> theta;
	vector<double> psi;
	vector<double> S1;
	vector<double> S2_num;
	vector<double> S2_den;
*/
	vector< vector<double> >	_data;
	vector<unsigned long int> number_density;

	double	angmax;
	double	angmin;
	double	angres;
	int		angbins;

	void OutputData ();
	//void OutputStatus () const;

};

#endif
