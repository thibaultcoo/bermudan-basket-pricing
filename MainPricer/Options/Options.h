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
		// Constructors
		Options();
		Options(StochasticProcess* process, double strike, std::vector<double> rates, double maturity);

		// Virtual functions for pricing using different methods
		virtual double price(int NbSim) = 0;
		virtual double priceAntithetic(int NbSim) = 0;
		virtual double priceControlVariate(int NbSim) = 0;
		virtual double priceQuasiMC(int NbSim) = 0;

		// Functions for statistical analysis
		double variance();
		double mean();
		std::vector<double> confidenceInterval(double alpha = 0.99);
};