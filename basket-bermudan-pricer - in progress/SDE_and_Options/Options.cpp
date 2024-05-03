#include "pch.h"
#include < algorithm >
#include <numeric>
#include "Options.h"

// Default Constructor
Options::Options() { }

// Initializer Constructor
Options::Options(StochasticProcess* process, double strike, std::vector<double> rates, double maturity)
	: process(process), strike(strike), rates(rates), maturity(maturity), paths_prices(1) { }

// Compute the mean of the simulated prices
double Options::mean()
{
	double nbSim = paths_prices.size();
	double sum = 0;
	for (size_t i = 0; i < nbSim; i++)
		sum += paths_prices[i];

	return sum / nbSim;
}

// Compute the variance of the simulated prices
double Options::variance()
{
	double nbSim = paths_prices.size();
	double mean = Options::mean();
	double sum = 0;
	for (size_t i = 0; i < nbSim; i++)
		sum += pow((paths_prices[i] - mean), 2);

	return sum / nbSim;
}

// Compute the Confidence Interval for the mean of the simulated prices
std::vector<double> Options::confidenceInterval(double alpha)
{
	if (alpha != 0.95 && alpha != 0.99)
		throw std::invalid_argument("Alpha should be 0.95 or 0.99");

	double z = (alpha == 0.99) ? 2.576 : 1.96;
	double mean_value = mean();
	double variance_value = variance();
	double stddev = std::sqrt(variance_value / paths_prices.size());

	std::vector<double> CI = { mean_value - z * stddev, mean_value + z * stddev };
	return CI;
}