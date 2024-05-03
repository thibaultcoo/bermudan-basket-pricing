#include "pch.h"
#include "EuropeanBasket.h"
#include "Tools.h"
#include <numeric>


EuropeanBasket::EuropeanBasket(StochasticProcess* process, const Instance& Inst)
	: Options(process, Inst.getStrike(), Inst.getRates(), Inst.getMatu()), weights(Inst.getWeights()), spots(Inst.getSpots()), corrMatrix(Inst.getVariances()) { }

// Pricing Methods
double EuropeanBasket::price(int nbSim)
{
	double sum = 0.;

	paths_prices.clear();
	paths_prices.resize(nbSim);

	for (int n = 0; n < nbSim; ++n)
	{
		process->Simulate(0, maturity, maturity * 365);
		std::vector<double> final_spots = process->Get_ValueND(maturity);
		double weighted_spots = std::inner_product(std::begin(weights), std::end(weights), std::begin(final_spots), 0.0);
		double path_price = std::max(weighted_spots - strike, 0.);
		sum += path_price;
		paths_prices[n] = path_price;
	}

	double price = std::exp(-rates[0] * maturity) * (sum / nbSim);

	return price;
}

double EuropeanBasket::priceAntithetic(int nbSim)
{
	double sum = 0.;

	paths_prices.clear();
	paths_prices.resize(nbSim);

	for (int n = 0; n < nbSim; ++n)
	{
		process->Simulate_Antithetic(0, maturity, maturity * 365);
		std::vector<double> final_spots = process->Get_ValueND(maturity);
		std::vector<double> final_spots_anti = process->Get_ValueND_antithetic(maturity);
		double weighted_spots = std::inner_product(std::begin(weights), std::end(weights), std::begin(final_spots), 0.0);
		double weighted_spots_anti = std::inner_product(std::begin(weights), std::end(weights), std::begin(final_spots_anti), 0.0);
		double path_price = std::max(weighted_spots - strike, 0.);
		double path_price_anti = std::max(weighted_spots_anti - strike, 0.);
		sum += (path_price + path_price_anti);
		paths_prices[n] = (path_price + path_price_anti) / 2;
	}

	double price = std::exp(-rates[0] * maturity) * (sum / (nbSim * 2));

	return price;
}

double EuropeanBasket::priceControlVariate(int nbSim)
{
	double sum = 0.;

	paths_prices.clear();
	paths_prices.resize(nbSim);

	double control_variate_payoff = compute_expected_value_control_variate(spots, weights, strike, rates[0], corrMatrix, maturity);

	for (int n = 0; n < nbSim; ++n)
	{
		process->Simulate(0, maturity, maturity * 365);
		std::vector<double> final_spots = process->Get_ValueND(maturity);
		std::vector<double> log_spots(spots.size());
		for (int j = 0; j < spots.size(); j++)
			log_spots[j] = log(final_spots[j]);

		double weighted_spots = std::inner_product(std::begin(weights), std::end(weights), std::begin(final_spots), 0.0);
		double weighted_log_spots = std::inner_product(std::begin(weights), std::end(weights), std::begin(log_spots), 0.0);

		double payoff = std::max(weighted_spots - strike, 0.);
		double log_payoff = std::max(exp(weighted_log_spots) - strike, 0.);
		double control_adjusted_payoff = payoff - log_payoff + control_variate_payoff;

		sum += control_adjusted_payoff;
		paths_prices[n] = control_adjusted_payoff;
	}

	double price = std::exp(-rates[0] * maturity) * (sum / nbSim);

	return price;
}

double EuropeanBasket::priceQuasiMC(int nbSim)
{
	double sum = 0.;

	paths_prices.clear();
	paths_prices.resize(nbSim);

	for (int n = 0; n < nbSim; ++n)
	{
		process->Simulate_VDC(0, maturity, maturity * 365, n, nbSim);
		std::vector<double> final_spots = process->Get_ValueND(maturity);
		double weighted_spots = std::inner_product(std::begin(weights), std::end(weights), std::begin(final_spots), 0.0);
		double path_price = std::max(weighted_spots - strike, 0.);
		sum += path_price;
		paths_prices[n] = path_price;
	}

	double price = std::exp(-rates[0] * maturity) * (sum / nbSim);

	return price;
}