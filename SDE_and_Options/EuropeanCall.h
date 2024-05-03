#pragma once
#include "Option.h"

//Pricing european call with Monte-Carlo simulations

class EuropeanCall : public Option
{
public:
	EuropeanCall();
	EuropeanCall(StochasticProcess* _process, double _K, std::vector<double> _r, double _T);

	double ComputePrice(int NbSim, bool antithetic = false);
	double ComputePrice_ControlVariate(int NbSim);
	double ComputePrice_VDC(int NbSim);
};