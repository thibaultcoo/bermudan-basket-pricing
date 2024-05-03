#include "Options.h"

Options::Options() { }

Options::Options(StochasticProcess* process, double strike, std::vector<double> rates, double maturity)
	: process(process), strike(strike), rates(rates), maturity(maturity), paths_prices(1) { }

double Options::mean()
{
	double nbSim = paths_prices.size();
	double sum = 0;
	for (size_t i = 0; i < nbSim; i++)
		sum += paths_prices[i];

	return sum / nbSim;
}

double Options::variance()
{
	double nbSim = paths_prices.size();
	double mean = Options::mean();
	double sum = 0;
	for (size_t i = 0; i < nbSim; i++)
		sum += pow((paths_prices[i] - mean), 2);

	return sum / nbSim;
}

std::vector<double> Options::confidenceInterval(double alpha)
{
	std::vector<double> CI(2);
	double z = 0.;

	if (alpha == 0.99)
		z = 2.576;
	else if (alpha == 0.95)
		z = 1.96;
	else
		throw std::exception("Alpha should be 0.99 or 0.95");

	double lower_bound = Options::mean() - z * sqrt(Options::variance() / paths_prices.size());
	double upper_bound = Options::mean() + z * sqrt(Options::variance() / paths_prices.size());
	CI[0] = lower_bound;
	CI[1] = upper_bound;

	return CI;
}