#ifndef MORITA2002_H_
#define MORITA2002_H_

#include <iostream>
#include "../utility.h"
#include "../watersystem.h"
#include "../moritasfg2002.h"

template <class T>
class Morita2002 : public WaterSystem<T> {

public:
	Morita2002 (const WaterSystemParams& params, const std::string filepath, const VecR& size, const std::string wannierpath)
	:
		WaterSystem<T>(params)
	{
		this->sys = new T(filepath, size, wannierpath);
	}

	Morita2002 (const WaterSystemParams& params, const std::string prmtop, const std::string mdcrd, const std::string mdvel = "")
	:
		WaterSystem<T>(params)
	{
		this->sys = new T(prmtop, mdcrd, mdvel);
	}

	void Analyze ();

private:
	VecR	_M;		// system polarization
	MatR	_alpha;	// system polarizability (susceptibility?)

	MoritaSFG	_sfg;

};

template <class T>
void Morita2002<T>::Analyze () {

	for (this->timestep = 0; this->timestep < this->timesteps; this->timestep++) {
		//cout << "step = " << step << endl;

		// find all the water molecules in the system
		/*
		RUN (this->sys->Molecules()) {
			Molecule * mol = this->sys->Molecules(i);
			if (mol->Name() != "h2o") {
				//mol->Print();
				std::cout << mol->Name() << " ";
			}
		}
			std::cout << std::endl;
		*/
		this->FindWaters ();
		this->SliceWaters (45.0, 70.0);
		//cout << this->int_mols.size() << endl;

		// calculate the polarization at time t=0
		if (!this->timestep)
			_M = _sfg.CalcTotalPolarization(this->int_mols);
			//_M.Print();

		// calculate the polarizability of the system at time t
		_alpha = _sfg.CalcTotalPolarizability(this->int_mols);

		// output the data so far
		double valx = 0.0;
		double valy = 0.0;
		double valz = 0.0;

		valx += _alpha.Index(1,1);
		valx += _alpha.Index(2,2);
		valx /= 2.0;
		valx *= _M[x];

		valy += _alpha.Index(0,0);
		valy += _alpha.Index(2,2);
		valy /= 2.0;
		valy *= _M[y];

		valz += _alpha.Index(0,0);
		valz += _alpha.Index(1,1);
		valz /= 2.0;
		valz *= _M[z];

		printf ("%d\t%d\t% 20.8f% 20.8f% 20.8f\n",
			this->timestep, (int)this->int_mols.size(), valx, valy, valz);
		fflush(stdout);

		this->sys->LoadNext();
	}

return;
}


#endif
