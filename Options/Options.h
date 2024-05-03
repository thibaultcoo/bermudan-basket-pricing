#pragma once
#include "RandomGenerator.h"
#include "StochasticProcess.h"


class Options
{
	protected:
		StochasticProcess* process;
		double strike;
		double maturity;
		std::vector<double> rates;
		std::vector<double> paths_prices;

	public:
		Options();
		Options(StochasticProcess* process, double strike, std::vector<double> rates, double maturity);
		virtual double price(int NbSim) = 0;
		virtual double priceAntithetic(int NbSim) = 0;
		virtual double priceControlVariate(int NbSim) = 0;
		virtual double priceVDC(int NbSim) = 0;
		double variance();
		double mean();
		std::vector<double> confidenceInterval(double alpha=0.99);
};

double fp(double x, double p);
