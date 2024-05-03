#pragma once
#include "Option.h"
#include <eigen-3.4.0/Eigen/Dense>

//Bermudan Call pricing

class BermudanCall : public Option
{

public:
	BermudanCall();
	BermudanCall(StochasticProcess* _process, double _K, std::vector<double> _r, double _T, std::vector<double> _exeDates, int _L = 5);
	double ComputePrice(int NbSim, bool antithetic = false);
	double ComputePrice_ControlVariate(int NbSim);
	double ComputePrice_VDC(int NbSim);

private:
	std::vector<double> exeDates;
	int L;
};