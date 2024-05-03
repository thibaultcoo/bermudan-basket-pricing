#pragma once
#include "RandomGenerator.h"
#include "StochasticProcess.h"

//Option main classes, EuropeanCall, BermudanCall,etc.. will inherit

class Option
{
public:
	Option();
	Option(StochasticProcess* _process, double _K, std::vector<double> _r, double _T);
	virtual double ComputePrice(int NbSim, bool antithetic = false) = 0;
	virtual double ComputePrice_ControlVariate(int NbSim) = 0;
	virtual double ComputePrice_VDC(int NbSim) = 0;
	double calculate_variance();
	double calculate_mean();
	std::vector<double> calculate_ConfidenceInterval(double alpha = 0.99);

protected:
	double K;
	std::vector<double> r;
	double T;
	StochasticProcess* process;
	std::vector<double> v;
};

double fp(double x, double p);
